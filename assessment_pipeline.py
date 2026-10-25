"""One command: download sources -> process datasets -> validate saved outputs.

Run: .venv/Scripts/python.exe assessment_pipeline.py
Dependencies: requirements.txt. Use --use-case1 or --use-case2 to select a workflow.
Outputs: data/processed, results/pipeline, results/use_case1 and results/use_case2.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'data/raw'
PROCESSED = ROOT / 'data/processed'
OUTPUT = ROOT / 'results/pipeline'
EBI = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90502001-GCST90503000'
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153837'
ACCESSIONS = ('GCST90502762', 'GCST90502763')
CHUNK_SIZE = 200_000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fetch(url, target, expected_md5=None):
    """Reuse complete files, resume partial transfers, and verify published hashes."""
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        partial = target.with_name(target.name + '.part')
        for attempt in range(1, 7):
            offset = partial.stat().st_size if partial.exists() else 0
            print(f'Download {target.name}: {offset:,} bytes saved; attempt {attempt}', flush=True)
            request = Request(url, headers={'Range': f'bytes={offset}-', 'Accept-Encoding': 'identity'})
            try:
                with urlopen(request, timeout=60) as response:
                    if response.status == 206:
                        match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
                        require(match is not None and int(match[1]) == offset, 'Invalid resumed response')
                        total, mode = int(match[3]), 'ab'
                    else:
                        require(response.status == 200, 'Unexpected download response')
                        total = int(response.headers.get('Content-Length', 0))
                        offset, mode = 0, 'wb'
                    reported = offset
                    with partial.open(mode) as output:
                        while block := response.read(65536):
                            output.write(block)
                            offset += len(block)
                            if offset - reported >= 32 * 1024 * 1024:
                                print(f'  {offset:,} / {total:,} bytes', flush=True)
                                reported = offset
                    require(not total or offset == total, 'Truncated download')
                break
            except HTTPError as error:
                if error.code == 416 and expected_md5 and partial.exists():
                    with partial.open('rb') as source:
                        if hashlib.file_digest(source, 'md5').hexdigest() == expected_md5:
                            break
                raise
            except (OSError, ValueError) as error:
                if attempt == 6:
                    raise
                print(f'  Retrying: {error}', flush=True)
        partial.replace(target)
    with target.open('rb') as source:
        md5 = hashlib.file_digest(source, 'md5').hexdigest()
    if expected_md5:
        require(md5 == expected_md5, f'Checksum mismatch: {target.name}')
    return {'file': str(target.relative_to(ROOT)), 'url': url,
            'bytes': target.stat().st_size, 'md5': md5,
            'published_md5_verified': expected_md5 is not None}


def download_sources(include_expression=True):
    records = []
    for tail in (('matrix/GSE153837_series_matrix.txt.gz', 'suppl/GSE153837_Raw_Data.txt.gz') if include_expression else ()):
        records.append(fetch(f'{GEO}/{tail}', RAW / tail.rsplit('/', 1)[-1]))
    metadata = ROOT / 'data/metadata'
    for accession in ACCESSIONS:
        base = f'{EBI}/{accession}/harmonised'
        for name in (f'{accession}.h.tsv.gz-meta.yaml', f'{accession}.running.log'):
            records.append(fetch(f'{base}/{name}', metadata / name))
        checksums = metadata / f'{accession}_harmonised_md5sum.txt'
        records.append(fetch(f'{base}/md5sum.txt', checksums))
        name = f'{accession}.h.tsv.gz'
        checksum = next(line.split()[0] for line in checksums.read_text().splitlines()
                        if line.split()[-1].lstrip('*') == name)
        records.append(fetch(f'{base}/{name}', RAW / name, checksum))
    (OUTPUT / 'downloads.json').write_text(json.dumps(records, indent=2), encoding='utf-8')


def process_expression():
    """Keep deposited scale; filter probes using assay detection and variation."""
    with gzip.open(RAW / 'GSE153837_series_matrix.txt.gz', 'rt') as source:
        text = source.read()
    header, body = text.split('!series_matrix_table_begin')
    fields = list(csv.reader(io.StringIO(header), delimiter='\t'))
    titles = next(row[1:] for row in fields if row and row[0] == '!Sample_title')
    identifiers = next(row[1:] for row in fields if row and row[0] == '!Sample_geo_accession')
    statuses = next(row[1:] for row in fields if len(row) > 1 and row[0] == '!Sample_characteristics_ch1'
                    and row[1].startswith('glucose metabolism status:'))
    samples = pd.DataFrame({'title': titles}, index=pd.Index(identifiers, name='sample_id'))
    samples['group'] = samples.title.str.extract(r'^\d+_(NGT|PD)_', expand=False)
    samples['time'] = samples.title.map(lambda x: 'fasting' if '_Fasting_' in x else '2h' if '_2_hour_' in x else None)
    samples['deposited_status'] = [value.split(': ', 1)[1] for value in statuses]
    expected = {('NGT', 'fasting'): 6, ('NGT', '2h'): 6, ('PD', 'fasting'): 9, ('PD', '2h'): 7}
    require(samples.groupby(['group', 'time']).size().to_dict() == expected, 'Incorrect sample groups')
    require(samples[['group', 'time']].notna().all().all(), 'Unknown group/time label')
    require(samples.deposited_status.eq(samples.group.map({'NGT': 'Normal', 'PD': 'Pre-diabetic'})).all(),
            'Sample titles disagree with deposited characteristics')
    expression = pd.read_csv(io.StringIO(body.split('!series_matrix_table_end')[0].strip()), sep='\t', index_col=0)
    require(expression.columns.tolist() == identifiers, 'Expression/sample order mismatch')
    require(expression.index.is_unique and samples.index.is_unique, 'Duplicate probe/sample identifiers')
    require(np.isfinite(expression.to_numpy()).all(), 'Invalid expression values')
    raw = pd.read_csv(RAW / 'GSE153837_Raw_Data.txt.gz', sep='\t', index_col=0)
    require(raw.index.equals(expression.index), 'Raw/processed probe mismatch')
    detection = raw[[title + '.Detection Pval' for title in titles]].copy()
    detection.columns = identifiers
    require(np.isfinite(detection.to_numpy()).all() and ((detection >= 0) & (detection <= 1)).all().all(),
            'Invalid detection p-values')
    # Retain probes detected in >= half the samples of at least one of the four
    # groups, allowing condition-specific expression. Do not test detection p-values as DE.
    detected = pd.Series(False, index=expression.index)
    for _, group in samples.groupby(['group', 'time']):
        detected |= (detection[group.index] < 0.01).sum(axis=1) >= int(np.ceil(len(group) / 2))
    variable = expression.var(axis=1) > 0
    keep = detected & variable
    cleaned = expression.loc[keep]
    require(len(cleaned) > 0, 'No expression probes retained')
    cleaned.to_csv(PROCESSED / 'expression_filtered.csv.gz', compression='gzip')
    samples.to_csv(PROCESSED / 'samples.csv')
    pd.DataFrame({'detected_in_any_group': detected, 'variable': variable, 'retained': keep}).to_csv(
        OUTPUT / 'probe_filtering.csv')
    samples['expression_sd'] = expression.std()
    samples['fraction_detected'] = (detection < 0.01).mean().to_numpy()
    samples.to_csv(OUTPUT / 'sample_quality.csv')
    check = pd.read_csv(PROCESSED / 'expression_filtered.csv.gz', index_col=0)
    require(check.shape == cleaned.shape and check.columns.tolist() == identifiers, 'Saved expression shape/order mismatch')
    require(np.allclose(check.to_numpy(), cleaned.to_numpy(), rtol=1e-12, atol=1e-12), 'Saved expression values changed')
    return {'input_probes': len(expression), 'retained_probes': len(cleaned), 'samples': len(samples),
            'excluded_samples': 0, 'validation': 'PASS',
            'scale': 'Deposited scale retained; no repeated log transformation or normalization',
            'review_needed': 'Preprocessing/centering and unusual PD sample dispersion remain unresolved; not final DE-ready data.'}


def process_gwas(accession):
    """Stream QC in chunks; keep valid harmonised autosomal polymorphic variants."""
    columns = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele',
               'beta', 'standard_error', 'effect_allele_frequency', 'p_value', 'rsid', 'hm_code']
    counts, chromosome_counts, missing = Counter(), Counter(), Counter()
    target = PROCESSED / f'{accession}_clean.csv.gz'
    temporary = target.with_name(target.name + '.part')
    previous_position, previous_alleles = None, set()
    top = pd.DataFrame()
    with gzip.open(temporary, 'wt', newline='') as output:
        first = True
        for chunk in pd.read_csv(RAW / f'{accession}.h.tsv.gz', sep='\t', usecols=columns,
                                 chunksize=CHUNK_SIZE, low_memory=False):
            counts['input'] += len(chunk)
            missing.update(chunk.isna().sum().to_dict())
            for column in ['chromosome', 'base_pair_location', 'beta', 'standard_error',
                           'effect_allele_frequency', 'p_value', 'hm_code']:
                chunk[column] = pd.to_numeric(chunk[column], errors='coerce')
            for column in ['effect_allele', 'other_allele']:
                chunk[column] = chunk[column].astype('string').str.upper()
            rules = [
                ('failed_harmonisation', chunk.hm_code.isin([5, 6, 10, 11])),
                ('invalid_or_non_autosomal_coordinates', chunk.chromosome.between(1, 22) & chunk.chromosome.mod(1).eq(0)
                 & chunk.base_pair_location.gt(0) & chunk.base_pair_location.mod(1).eq(0)),
                ('invalid_alleles', chunk.effect_allele.str.fullmatch('[ACGT]+', na=False)
                 & chunk.other_allele.str.fullmatch('[ACGT]+', na=False) & chunk.effect_allele.ne(chunk.other_allele)),
                ('invalid_effect_or_se', np.isfinite(chunk.beta) & np.isfinite(chunk.standard_error) & chunk.standard_error.gt(0)),
                ('invalid_p_value', chunk.p_value.gt(0) & chunk.p_value.le(1)),
                ('invalid_or_monomorphic_frequency', chunk.effect_allele_frequency.gt(0) & chunk.effect_allele_frequency.lt(1)),
            ]
            keep = pd.Series(True, index=chunk.index)
            for reason, valid in rules:
                valid = valid.fillna(False)
                counts[reason] += int((keep & ~valid).sum())
                keep &= valid
            cleaned = chunk.loc[keep].copy()
            if cleaned.empty:
                continue
            cleaned[['chromosome', 'base_pair_location']] = cleaned[['chromosome', 'base_pair_location']].astype('int64')
            positions = pd.MultiIndex.from_frame(cleaned[['chromosome', 'base_pair_location']])
            require(positions.is_monotonic_increasing, 'GWAS input is not position sorted')
            require(previous_position is None or positions[0] >= previous_position, 'GWAS chunk order is invalid')
            # Sorted coordinates let us detect duplicates across chunk boundaries
            # while keeping only the last coordinate's allele pairs in memory.
            key_columns = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele']
            duplicate = cleaned.duplicated(key_columns)
            if previous_position is not None:
                at_boundary = (cleaned.chromosome == previous_position[0]) & (cleaned.base_pair_location == previous_position[1])
                duplicate |= at_boundary & pd.Series(list(zip(cleaned.effect_allele, cleaned.other_allele)), index=cleaned.index).isin(previous_alleles)
            counts['duplicate_variant'] += int(duplicate.sum())
            last_position = positions[-1]
            last_rows = cleaned.loc[(cleaned.chromosome == last_position[0]) & (cleaned.base_pair_location == last_position[1])]
            last_alleles = set(zip(last_rows.effect_allele, last_rows.other_allele))
            previous_alleles = previous_alleles | last_alleles if last_position == previous_position else last_alleles
            previous_position = last_position
            cleaned = cleaned.loc[~duplicate]
            cleaned['minor_allele_frequency'] = np.minimum(cleaned.effect_allele_frequency, 1 - cleaned.effect_allele_frequency)
            counts['retained'] += len(cleaned)
            counts['genome_wide_significant'] += int((cleaned.p_value < 5e-8).sum())
            chromosome_counts.update(cleaned.chromosome.value_counts().to_dict())
            cleaned.to_csv(output, index=False, header=first)
            first = False
            top = pd.concat([top, cleaned.nsmallest(100, 'p_value')]).nsmallest(100, 'p_value')
            print(f'{accession}: processed {counts["input"]:,}; retained {counts["retained"]:,}', flush=True)
    excluded = sum(counts[reason] for reason in [r[0] for r in rules]) + counts['duplicate_variant']
    require(counts['input'] == counts['retained'] + excluded, 'GWAS counts do not reconcile')
    require(counts['retained'] > 0, 'No GWAS variants retained')
    # Read the entire output back, validating gzip integrity and row counts.
    saved = sum(len(block) for block in pd.read_csv(temporary, compression='gzip', usecols=['p_value'], chunksize=CHUNK_SIZE))
    require(saved == counts['retained'], 'Saved GWAS row count mismatch')
    temporary.replace(target)
    top.to_csv(OUTPUT / f'{accession}_top100.csv', index=False)
    return {'counts': dict(counts), 'chromosome_counts': {str(k): v for k, v in chromosome_counts.items()},
            'input_missing_by_column': dict(missing), 'validation': 'PASS',
            'notes': 'No MAF cutoff beyond polymorphism; no imputation metrics supplied. Palindromic variants retained. Duplicate key uses oriented effect/other alleles. Top100 is statistical ranking, not independent loci.'}


def analyse_use_case1(qc):
    """Complete GWAS figures, descriptive loci, BMI sensitivity and scientific report."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import chi2, beta as beta_distribution

    out = ROOT / 'results/use_case1'
    out.mkdir(parents=True, exist_ok=True)
    cutoff = 5e-8
    key = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele']
    tops = {a: pd.read_csv(OUTPUT / f'{a}_top100.csv') for a in ACCESSIONS}
    union = pd.concat([d.loc[d.p_value < cutoff] for d in tops.values()]).drop_duplicates(key)
    require(not union.empty, 'No significant signals; report template needs review')
    # Distance grouping is descriptive; no LD independence or causal claims.
    loci = []
    for row in union.sort_values('p_value').itertuples():
        if not any(row.chromosome == c and abs(row.base_pair_location - p) <= 500_000 for c, p in loci):
            loci.append((row.chromosome, row.base_pair_location))
    union_keys = set(map(tuple, union[key].to_numpy()))
    summaries, selected = {}, {}
    for accession in ACCESSIONS:
        p_arrays, maf_arrays, scatter, regions, selected_rows, significant = [], [], [], [], [], []
        chromosome_max = Counter()
        rng = np.random.default_rng(2026)
        count = 0
        for block in pd.read_csv(PROCESSED / f'{accession}_clean.csv.gz', chunksize=CHUNK_SIZE):
            require(np.isfinite(block[['beta', 'standard_error', 'p_value', 'effect_allele_frequency']]).all().all(), 'Nonfinite GWAS values')
            require(block.p_value.between(0, 1, inclusive='right').all() and block.standard_error.gt(0).all(), 'Invalid association data')
            count += len(block)
            p_arrays.append(block.p_value.to_numpy())
            maf_arrays.append(block.minor_allele_frequency.to_numpy())
            keep_plot = (rng.random(len(block)) < 0.02) | block.p_value.lt(1e-5)
            scatter.append(block.loc[keep_plot, ['chromosome', 'base_pair_location', 'p_value']])
            for chromosome, position in loci:
                region = block.loc[block.chromosome.eq(chromosome) & block.base_pair_location.sub(position).abs().le(500_000)]
                if len(region):
                    regions.append(region)
            members = pd.Series(list(map(tuple, block[key].to_numpy())), index=block.index).isin(union_keys)
            selected_rows.append(block.loc[members])
            significant.append(block.loc[block.p_value < cutoff])
            for chromosome, maximum in block.groupby('chromosome').base_pair_location.max().items():
                chromosome_max[int(chromosome)] = max(chromosome_max[int(chromosome)], int(maximum))
        require(count == qc[accession]['counts']['retained'], 'Clean GWAS count differs from validated processing run')
        pvalues, maf = np.concatenate(p_arrays), np.concatenate(maf_arrays)
        common = maf >= 0.01
        lam = float(chi2.isf(np.median(pvalues), 1) / chi2.ppf(0.5, 1))
        lam_common = float(chi2.isf(np.median(pvalues[common]), 1) / chi2.ppf(0.5, 1))
        hits = pd.concat(significant, ignore_index=True)
        selected[accession] = pd.concat(selected_rows, ignore_index=True)
        for table in (hits, selected[accession]):
            table['odds_ratio'] = np.exp(table.beta)
            table['ci_lower'] = np.exp(table.beta - 1.95996398454 * table.standard_error)
            table['ci_upper'] = np.exp(table.beta + 1.95996398454 * table.standard_error)
        hits.to_csv(out / f'{accession}_significant_variants.csv', index=False)
        summaries[accession] = {'variants': count, 'significant_variants': len(hits),
            'lambda_gc': lam, 'lambda_gc_maf_ge_0_01': lam_common,
            'maf_below_0_01': int((~common).sum()),
            'significant_maf_ge_0_01': int(((pvalues < cutoff) & common).sum()),
            'min_p': float(pvalues.min()), 'validation': 'PASS'}
        # QQ statistics use ALL retained variants. Only plotted ranks are reduced.
        ordered = np.sort(pvalues)
        ranks = np.unique(np.r_[np.arange(1, min(1001, count + 1)),
                                   np.geomspace(1, count, 6000).astype(int)])
        expected = -np.log10((ranks - 0.5) / count)
        observed = -np.log10(ordered[ranks - 1])
        fig, ax = plt.subplots(figsize=(6, 6))
        lo = -np.log10(beta_distribution.ppf(.025, ranks, count + 1 - ranks))
        hi = -np.log10(beta_distribution.ppf(.975, ranks, count + 1 - ranks))
        order = np.argsort(expected)
        ax.fill_between(expected[order], hi[order], lo[order], color='lightgray', label='Pointwise 95% null envelope')
        ax.scatter(expected, observed, s=5, color='#176b99')
        end = max(expected.max(), observed.max()) + .2
        ax.plot([0, end], [0, end], '--', color='black', linewidth=1)
        ax.set(xlabel='Expected -log10(p)', ylabel='Observed -log10(p)',
               title=f'{accession}: Q–Q; λGC={lam:.3f}')
        ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(out / f'{accession}_qq.png', dpi=180); plt.close(fig)
        points = pd.concat(scatter, ignore_index=True)
        offsets, ticks, offset = {}, [], 0
        for chromosome in sorted(chromosome_max):
            offsets[chromosome] = offset
            ticks.append((offset + chromosome_max[chromosome] / 2, str(chromosome)))
            offset += chromosome_max[chromosome] + 5_000_000
        fig, ax = plt.subplots(figsize=(13, 4.5))
        for chromosome, group in points.groupby('chromosome'):
            ax.scatter(group.base_pair_location + offsets[chromosome], -np.log10(group.p_value),
                       s=3, color=('#176b99' if chromosome % 2 else '#8899aa'), rasterized=True)
        ax.axhline(-np.log10(cutoff), color='firebrick', linestyle='--', linewidth=1)
        for chromosome, position in loci:
            lead = hits.loc[hits.chromosome.eq(chromosome)].nsmallest(1, 'p_value')
            if len(lead):
                r = lead.iloc[0]
                ax.annotate(r.rsid, (r.base_pair_location + offsets[chromosome], -np.log10(r.p_value)),
                            xytext=(0, 6), textcoords='offset points', fontsize=8, ha='center')
        ax.set(xticks=[t[0] for t in ticks], xticklabels=[t[1] for t in ticks],
               xlabel='Chromosome (GRCh38)', ylabel='-log10(p)', title=f'{accession}: Manhattan plot')
        ax.set_ylim(0, max(9.5, -np.log10(pvalues.min()) + 1))
        fig.tight_layout(); fig.savefig(out / f'{accession}_manhattan.png', dpi=180); plt.close(fig)
        regional = pd.concat(regions, ignore_index=True).drop_duplicates(key)
        regional.to_csv(out / f'{accession}_regional_variants.csv.gz', index=False, compression='gzip')
        fig, axes = plt.subplots(1, len(loci), figsize=(6 * len(loci), 4), squeeze=False)
        for ax, (chromosome, position) in zip(axes[0], loci):
            r = regional.loc[regional.chromosome.eq(chromosome) & regional.base_pair_location.sub(position).abs().le(500_000)]
            ax.scatter(r.base_pair_location / 1e6, -np.log10(r.p_value), s=7, color='#176b99')
            ax.axhline(-np.log10(cutoff), color='firebrick', linestyle='--')
            ax.set(xlabel=f'Chromosome {chromosome} position (Mb, GRCh38)', ylabel='-log10(p)',
                   title='Regional associations; LD not estimated')
        fig.tight_layout(); fig.savefig(out / f'{accession}_regional.png', dpi=180); plt.close(fig)
        print(f'{accession}: {count:,} variants; {len(hits)} significant; lambda={lam:.4f}', flush=True)
    comparison = selected[ACCESSIONS[0]].merge(selected[ACCESSIONS[1]], on=key, suffixes=('_primary', '_bmi'), validate='one_to_one')
    require(len(comparison) == len(union), 'Significant union not fully allele-aligned between models')
    comparison['same_effect_direction'] = np.sign(comparison.beta_primary) == np.sign(comparison.beta_bmi)
    comparison['candidate_gene'] = comparison.chromosome.map({7: 'GCK', 10: 'HK1'}).fillna('Requires review')
    comparison['gene_evidence'] = 'Literature-supported locus annotation; not causal variant-to-gene proof'
    comparison.to_csv(out / 'bmi_sensitivity_comparison.csv', index=False)
    leads = []
    for chromosome, position in loci:
        candidate = comparison.loc[comparison.chromosome.eq(chromosome) & comparison.base_pair_location.sub(position).abs().le(500_000)]
        lead = candidate.nsmallest(1, 'p_value_primary').copy()
        lead['significant_variants_primary'] = int(candidate.p_value_primary.lt(cutoff).sum())
        lead['significant_variants_bmi'] = int(candidate.p_value_bmi.lt(cutoff).sum())
        leads.append(lead)
    leads = pd.concat(leads, ignore_index=True)
    leads.to_csv(out / 'prioritized_loci.csv', index=False)
    fig, ax = plt.subplots(figsize=(8, 4))
    for index, row in leads.iterrows():
        for delta, suffix, label, color in [(-.10, 'primary', 'Primary', '#176b99'), (.10, 'bmi', 'BMI-adjusted', '#d87922')]:
            estimate = row[f'odds_ratio_{suffix}']
            ax.errorbar(estimate, index + delta, xerr=[[estimate - row[f'ci_lower_{suffix}']], [row[f'ci_upper_{suffix}'] - estimate]],
                        fmt='o', color=color, label=label if index == 0 else None, capsize=3)
    ax.axvline(1, color='gray', linestyle='--')
    ax.set(yticks=range(len(leads)), yticklabels=[f'{r.rsid_primary} ({r.candidate_gene}), {r.effect_allele}' for r in leads.itertuples()],
           xlabel='Odds ratio per deposited effect allele (95% Wald CI)', title='Lead variant effects with and without BMI adjustment')
    ax.legend(); fig.tight_layout(); fig.savefig(out / 'lead_effects.png', dpi=180); plt.close(fig)
    summary = {'models': summaries, 'descriptive_locus_count': len(leads),
        'significant_union_count': len(union), 'all_aligned_effect_directions_agree': bool(comparison.same_effect_direction.all()),
        'genome_wide_threshold': cutoff, 'grouping_distance_bp': 500000,
        'plotting': 'Manhattan: seeded 2% background plus all p<1e-5; QQ: all variants used, reduced ranks displayed; regional: all variants within windows.',
        'status': 'PASS'}
    require(all(comparison.same_effect_direction), 'Opposite effects require report review')
    (out / 'analysis_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_use_case1_report(summary, leads, comparison, qc, out)
    expected_files = ['use_case1_report.md', 'use_case1_report.html', 'prioritized_loci.csv', 'bmi_sensitivity_comparison.csv', 'lead_effects.png']
    expected_files += [f'{a}_{suffix}' for a in ACCESSIONS for suffix in ['manhattan.png', 'qq.png', 'regional.png', 'significant_variants.csv']]
    require(all((out / name).is_file() and (out / name).stat().st_size > 0 for name in expected_files), 'Missing use case 1 outputs')
    (out / 'validation.json').write_text(json.dumps({'status': 'PASS', 'checked_outputs': expected_files,
        'checks': ['raw GWAS MD5 (download step)', 'QC counts reconcile', 'complete output row counts',
                   'all p-values used for QQ', 'significant union allele-aligned', 'effect directions agree', 'required output files exist'],
        'completed_utc': datetime.now(timezone.utc).isoformat()}, indent=2), encoding='utf-8')


def write_use_case1_report(summary, leads, comparison, qc, out):
    def effect(row, suffix):
        return f'{row[f"odds_ratio_{suffix}"]:.3f} ({row[f"ci_lower_{suffix}"]:.3f}–{row[f"ci_upper_{suffix}"]:.3f})'
    table = '| Lead variant | GRCh38 | Effect/other | Candidate gene | Primary OR (95% CI) | Primary p | BMI OR (95% CI) | BMI p |\n|---|---|---|---|---|---|---|---|\n'
    for _, r in leads.iterrows():
        table += f'| {r.rsid_primary} | {r.chromosome}:{r.base_pair_location} | {r.effect_allele}/{r.other_allele} | {r.candidate_gene} | {effect(r,"primary")} | {r.p_value_primary:.3g} | {effect(r,"bmi")} | {r.p_value_bmi:.3g} |\n'
    metrics = '| Dataset | Retained variants | p < 5×10⁻⁸ | λGC | MAF < 1% |\n|---|---:|---:|---:|---:|\n'
    for accession, m in summary['models'].items():
        metrics += f'| {accession} | {m["variants"]:,} | {m["significant_variants"]} | {m["lambda_gc"]:.4f} | {m["maf_below_0_01"]:,} |\n'
    text = f'''# Use Case 1 — GWAS summary-statistics analysis of prediabetes

## Objective

Identify prediabetes-associated variants, prioritize reproducible signals across BMI adjustment, and assess biological relevance. Prediabetes is relevant to early metabolic dysfunction; these case–control statistics concern prevalent prediabetes, not longitudinal progression to diabetes.

## Methodology

### Dataset exploration (Task 1)

The supplied accessions describe Qatari whole-genome sequencing on GRCh38. GCST90502762 adjusts for age, sex and PC1–PC10; GCST90502763 additionally adjusts for BMI. Both represent the same reported cohort, so BMI adjustment is a sensitivity comparison, not replication. [Catalog metadata](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90502001-GCST90503000/GCST90502762/GCST90502762.tsv-meta.yaml)

The publication specifies 2,267 cases and 8,665 controls. Prediabetes was defined using HbA1c 5.7–6.4% or fasting glucose 5.6–6.9 mmol/L in participants without declared diabetes/treatment. The original association model used SAIGE. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12232356/)

The harmonised files contain chromosome, base-pair position, effect/other alleles, beta, standard error, effect-allele frequency, p-value, rs identifiers, harmonisation flags and variant identifiers. The processing output retains the analysis columns and adds minor-allele frequency. Effects are log odds for the deposited effect allele; beta/standard-error Wald intervals are reported alongside the original association p-values. SAIGE p-values need not equal simple Wald p-values.

### Quality assessment (Task 2)

Downloads passed published MD5 checksums. Numeric fields, positive integral autosomal coordinates, distinct DNA alleles, finite effects, positive standard errors, 0<p≤1 and 0<EAF<1 were checked. Only successful harmonisation codes 5, 6, 10 and 11 were accepted. Exact duplicates were checked using chromosome, position and oriented alleles, including across chunk boundaries. All retained output rows were read back; counts reconcile. No additional local exclusions were needed in either supplied harmonised file.

Deposited harmonisation logs report 149,810 failed sites for GCST90502762 and 149,781 for GCST90502763; these failed sites are already absent from the downloaded harmonised files. Those upstream exclusions are distinct from local QC. [Primary harmonisation log](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90502001-GCST90503000/GCST90502762/harmonised/GCST90502762.running.log), [BMI harmonisation log](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90502001-GCST90503000/GCST90502763/harmonised/GCST90502763.running.log)

Polymorphic variants were retained without another MAF cutoff; the table below audits low-frequency content. The publication reports upstream MAF, call-rate and Hardy–Weinberg filtering. Individual-level QC cannot be repeated from summary statistics. Palindromic variants were retained for within-study analysis; no cross-strand guessing was used. Genomic inflation was calculated from all p-values as χ²₁⁻¹(1−median(p))/median(χ²₁).

### Statistical prioritization (Task 3)

Genome-wide significance was predefined as p<5×10⁻⁸ for the primary model. BMI-adjusted results assess sensitivity. Variants are ranked by original p-value; signals within 500 kb of a lead on the same chromosome form descriptive groups. These groups do not establish LD independence. Odds ratios are exp(beta), with 95% Wald intervals exp(beta±1.96×SE). The union of significant variants was matched by position and both oriented alleles between models. No meta-analysis was performed because participants overlap.

## Results

{metrics}

There are {summary['significant_union_count']} distinct significant variants across the two models, forming {summary['descriptive_locus_count']} separated genomic regions. All aligned variants have consistent effect directions. BMI adjustment strengthens the lead p-values, but smaller p-values alone do not prove a BMI-independent causal mechanism. The four primary significant variants comprise three on chromosome 7 and one on chromosome 10; the BMI model adds three nearby variants rather than an additional region.

{table}

The GCK-region A allele at rs2908293 increases prediabetes odds. The HK1-region G allele at rs72805692 decreases odds; its alternative A allele is therefore risk-increasing. The published study uses the opposite allele orientation at that HK1 lead. This is an orientation difference, not a conflicting result. The present signals agree with the study's reported GCK/HK1 regions; they are a reanalysis, not a new discovery. [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12232356/)

### Figures and interpretation (Task 4)

1. **Manhattan plots:** genomic locations versus −log10(p), with the red line at 5×10⁻⁸. Two regions exceed this threshold. A seeded 2% background plus every p<10⁻⁵ point is displayed; all variants contribute to testing and counts.
2. **Q–Q plots:** observed versus null p-value distributions, with a pointwise 95% envelope. λGC summarizes bulk calibration, while the upper tail shows the strongest associations. The envelope assumes independent tests and is approximate under LD. Calibration near one cannot rule out all confounding.
3. **Regional plots:** every retained variant within ±500 kb of each lead. Nearby peaks indicate local association structure; the plots contain no LD colouring and do not constitute fine-mapping.
4. **Effect plot:** lead odds ratios and confidence intervals for both models, using deposited effect alleles. Consistency supports robustness to model adjustment; the two estimates are correlated.

![Primary Manhattan](GCST90502762_manhattan.png)
![Primary QQ](GCST90502762_qq.png)
![BMI Manhattan](GCST90502763_manhattan.png)
![BMI QQ](GCST90502763_qq.png)
![Primary regional associations](GCST90502762_regional.png)
![BMI regional associations](GCST90502763_regional.png)
![Lead effects](lead_effects.png)

## Biological and clinical interpretation (Task 5)

**GCK:** Glucokinase phosphorylates glucose and supports pancreatic glucose-stimulated insulin secretion and hepatic glucose handling. These functions make the locus biologically plausible for a prediabetes phenotype. Rare pathogenic GCK mutations can cause monogenic disease, but that does not make these common associated variants pathogenic or diagnostic. The lead coordinate lies inside the current GRCh38 GCK gene interval; a gene-body overlap does not identify the causal transcript or regulatory mechanism. [NCBI GCK](https://www.ncbi.nlm.nih.gov/gene/2645)

**HK1:** Hexokinase 1 connects glucose phosphorylation with erythrocyte biology. HK1-region associations with HbA1c can reflect blood-cell mechanisms as well as glycaemia. Consequently, an association with an HbA1c-based case definition should not automatically be interpreted as altered insulin resistance. Previous research links HK1 variation to HbA1c and a proanemic state; that work concerns locus biology and does not prove the same mechanism for this lead variant. [HK1 primary study](https://diabetesjournals.org/diabetes/article/58/11/2687/15337/Genetic-Variant-in-HK1-Is-Associated-With-a), [HbA1c glycaemic/nonglycaemic mechanisms](https://pmc.ncbi.nlm.nih.gov/articles/PMC2992787/)

Candidate-gene assignments are literature-supported locus annotations. They are not evidence of variant pathogenicity, colocalisation or causal gene assignment. No pathway-enrichment result is claimed from only two candidate genes. The shared glucose-phosphorylation theme is a biological interpretation. These associations alone do not justify clinical screening changes, drug selection or individual risk prediction.

## Proposed next steps (Task 6)

1. **Independent replication:** test the lead alleles in another Qatari/Middle Eastern cohort, then other ancestries, using harmonised case definitions. This tests reproducibility and transferability.
2. **Separate phenotype components:** compare fasting-glucose-defined and HbA1c-defined prediabetes; examine erythrocyte indices and iron status, especially at HK1. This distinguishes glycaemic effects from measurement-related mechanisms.
3. **LD and fine-mapping:** use ancestry-matched genotypes or an appropriate LD reference for clumping, conditional association and credible sets. This resolves whether nearby significant variants tag one signal.
4. **Functional evidence:** annotate variants with Ensembl VEP and regulatory resources; test colocalisation with pancreatic-islet/liver eQTLs for GCK and erythroid evidence for HK1. Follow supported regulatory hypotheses with reporter or CRISPR experiments.
5. **Robustness and prospective validation:** evaluate age sensitivity, relatedness, sex and baseline BMI when individual data are available. Test progression in a separate longitudinal cohort; only then assess whether a validated score improves calibration, discrimination and clinical decisions beyond standard predictors.

## Limitations

Only summary statistics are available. Individual call rate, Hardy–Weinberg equilibrium, imputation quality, participant covariates and disease classification cannot be rechecked. Two analyses of the same participants do not provide replication. Distance grouping is not LD clumping; no ancestry-matched LD, conditional testing, credible sets or functional experiments were performed. Gene assignment and HK1 mechanism remain hypotheses. A fixed genome-wide threshold is conventional rather than a dataset-specific independent-test correction. Effect estimates are subject to winner's curse; modest sample size limits power. Ancestry and phenotype definition limit generalisability. This analysis cannot establish causality or predict progression from case–control statistics.

## Conclusion

The supplied prediabetes summary statistics identify two prioritized regions near GCK and HK1, with consistent allele-aligned effects after BMI adjustment. The results support further investigation of glucose regulation and HbA1c-related erythrocyte biology. Independent replication and functional evidence are required before clinical translation.

## Reproducibility

Code: [single Python pipeline](../../assessment_pipeline.py). Dependencies: [requirements](../../requirements.txt).

Fresh run: `python assessment_pipeline.py --use-case1`.
Reuse previously validated processed data: `python assessment_pipeline.py --use-case1 --reuse-processed`.

The pipeline downloads and verifies sources, processes data, generates all statistical tables/figures and this report, then writes [validation](validation.json). [Analysis metrics](analysis_summary.json), [lead loci](prioritized_loci.csv), [all significant model comparisons](bmi_sensitivity_comparison.csv).

The code repository link is supplied separately with the assessment submission.
'''
    (out / 'use_case1_report.md').write_text(text, encoding='utf-8')
    import base64
    import markdown
    html_body = markdown.markdown(text, extensions=['tables', 'fenced_code'])
    for picture in out.glob('*.png'):
        encoded = base64.b64encode(picture.read_bytes()).decode('ascii')
        html_body = html_body.replace(f'src="{picture.name}"', f'src="data:image/png;base64,{encoded}"')
    document = '''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Use Case 1 — Prediabetes GWAS</title>
<style>body{max-width:1100px;margin:40px auto;padding:0 24px;font:16px/1.6 system-ui;color:#203040}
h1,h2,h3{line-height:1.25;color:#15486b}h2{margin-top:2em}table{border-collapse:collapse;font-size:13px;width:100%}
th,td{border:1px solid #ccd7df;padding:8px;text-align:left}th{background:#edf4f8}img{max-width:100%;height:auto}
a{color:#176b99}code{background:#f1f4f6;padding:2px 5px} @media print{body{margin:0;font-size:11px}img,table{break-inside:avoid}h2,h3{break-after:avoid}}</style><body>'''
    (out / 'use_case1_report.html').write_text(document + html_body + '</body></html>', encoding='utf-8')


def run_use_case1(qc):
    out = ROOT / 'results/use_case1'
    out.mkdir(parents=True, exist_ok=True)
    validation = out / 'validation.json'
    validation.write_text(json.dumps({'status': 'RUNNING'}), encoding='utf-8')
    try:
        analyse_use_case1(qc)
    except Exception as error:
        validation.write_text(json.dumps({'status': 'FAIL', 'error': str(error)}, indent=2), encoding='utf-8')
        raise




def expression_annotations(archive_path):
    """Read the deposited R2 manifest in memory; never extract untrusted paths."""
    import tarfile
    with tarfile.open(archive_path) as archive:
        member = next(m for m in archive.getmembers() if '_R2_' in m.name and m.name.endswith('.txt.gz'))
        lines = gzip.decompress(archive.extractfile(member).read()).decode('utf-8', errors='replace').splitlines()
    section = '\n'.join(lines[lines.index('[Probes]') + 1:]).split('[Controls]')[0]
    annotation = pd.read_csv(io.StringIO(section), sep='\t', dtype=str)
    annotation = annotation.loc[annotation.Species.eq('Homo sapiens'), ['Probe_Id', 'Entrez_Gene_ID', 'Symbol', 'Definition']].drop_duplicates()
    annotation = annotation.loc[annotation.Entrez_Gene_ID.str.fullmatch(r'\d+', na=False) & annotation.Symbol.notna()]
    ambiguous = annotation.groupby('Probe_Id').Entrez_Gene_ID.nunique()
    annotation = annotation.loc[annotation.Probe_Id.isin(ambiguous.index[ambiguous.eq(1)])].drop_duplicates('Probe_Id')
    return annotation.set_index('Probe_Id')


def welch_expression(expression, samples, time):
    """Independent PD-minus-NGT test at one time; no equal-variance assumption."""
    from scipy.stats import ttest_ind, t
    from statsmodels.stats.multitest import multipletests
    pd_ids = samples.index[samples.group.eq('PD') & samples.time.eq(time)]
    ngt_ids = samples.index[samples.group.eq('NGT') & samples.time.eq(time)]
    require(len(pd_ids) >= 3 and len(ngt_ids) >= 3, 'Insufficient biological replicates')
    a, b = expression[pd_ids].to_numpy(), expression[ngt_ids].to_numpy()
    result = ttest_ind(a, b, axis=1, equal_var=False)
    difference = a.mean(axis=1) - b.mean(axis=1)
    standard_error = np.sqrt(a.var(axis=1, ddof=1) / len(pd_ids) + b.var(axis=1, ddof=1) / len(ngt_ids))
    valid = np.isfinite(result.pvalue) & np.isfinite(result.statistic) & (standard_error > 0)
    table = pd.DataFrame({'mean_PD': a.mean(axis=1), 'mean_NGT': b.mean(axis=1), 'mean_difference': difference,
                          'standard_error': standard_error, 't_statistic': result.statistic,
                          'welch_df': result.df, 'p_value': np.where(valid, result.pvalue, 1.0),
                          'testable': valid}, index=expression.index)
    # Retain untestable genes with p=1 so the tested family is never silently reduced.
    table['fdr_bh'] = multipletests(table.p_value, method='fdr_bh')[1]
    margin = t.ppf(.975, result.df) * standard_error
    table['ci_lower'], table['ci_upper'] = difference - margin, difference + margin
    table['significant'] = table.fdr_bh < .05
    table['direction'] = np.where(table.mean_difference > 0, 'higher_in_PD', 'lower_in_PD')
    table['n_PD'], table['n_NGT'] = len(pd_ids), len(ngt_ids)
    require(table.p_value.between(0, 1).all() and table.fdr_bh.between(0, 1).all(), 'Invalid DE probabilities')
    return table.sort_values(['fdr_bh', 'p_value'])


def quantile_normalize(matrix):
    """Quantile normalization with averaged ranks for ties."""
    from scipy.stats import rankdata
    values = matrix.to_numpy()
    reference = np.sort(values, axis=0).mean(axis=1)
    result = np.empty_like(values)
    for column in range(values.shape[1]):
        ranks = rankdata(values[:, column], method='average') - 1
        result[:, column] = np.interp(ranks, np.arange(len(reference)), reference)
    return pd.DataFrame(result, index=matrix.index, columns=matrix.columns)


def hallmark_scores(expression, annotation, samples, time, gene_sets):
    """Self-contained mean-z module tests; preserve samples as independent units."""
    from itertools import combinations
    from statsmodels.stats.multitest import multipletests
    ids = samples.index[samples.time.eq(time)]
    x = expression[ids]
    sd = x.std(axis=1)
    z = x.sub(x.mean(axis=1), axis=0).div(sd.replace(0, np.nan), axis=0)
    scores, overlaps = {}, {}
    symbols = annotation.Symbol
    for term, members in gene_sets.items():
        included = expression.index[symbols.isin(members) & sd.gt(0)]
        if len(included) >= 10:
            scores[term] = z.loc[included].mean(axis=0)
            overlaps[term] = sorted(set(symbols.loc[included]))
    score_frame = pd.DataFrame(scores, index=ids).T
    require(len(score_frame) > 0, 'No gene sets have adequate measured coverage')
    stats = welch_expression(score_frame, samples.loc[ids], time)
    # Enumerate every group allocation as an exchangeability sensitivity check.
    # Welch p-values remain primary because these samples show unequal dispersion.
    values = score_frame.to_numpy()
    n = len(ids)
    n_pd = int(samples.loc[ids].group.eq('PD').sum())
    masks = np.zeros((len(list(combinations(range(n), n_pd))), n))
    for index, selected in enumerate(combinations(range(n), n_pd)):
        masks[index, list(selected)] = 1
    sums_a = values @ masks.T
    sums_b = values.sum(axis=1, keepdims=True) - sums_a
    squares_a = (values ** 2) @ masks.T
    squares_b = (values ** 2).sum(axis=1, keepdims=True) - squares_a
    variance_a = np.maximum((squares_a - sums_a ** 2 / n_pd) / (n_pd - 1), 0)
    variance_b = np.maximum((squares_b - sums_b ** 2 / (n - n_pd)) / (n - n_pd - 1), 0)
    perm_t = (sums_a / n_pd - sums_b / (n - n_pd)) / np.sqrt(variance_a / n_pd + variance_b / (n - n_pd))
    observed = stats.loc[score_frame.index, 't_statistic'].to_numpy()
    exact = (np.abs(perm_t) >= np.abs(observed[:, None]) - 1e-10).mean(axis=1)
    stats['permutation_p_sensitivity'] = pd.Series(exact, index=score_frame.index)
    stats['permutation_fdr_sensitivity'] = multipletests(stats.permutation_p_sensitivity, method='fdr_bh')[1]
    stats['measured_genes'] = [len(overlaps[term]) for term in stats.index]
    stats['library_genes'] = [len(gene_sets[term]) for term in stats.index]
    stats['coverage'] = stats.measured_genes / stats.library_genes
    stats['measured_members'] = [';'.join(overlaps[term]) for term in stats.index]
    stats['permutations'] = len(masks)
    stats.index.name = 'hallmark'
    return stats, score_frame


def analyse_use_case2(out):
    import tarfile
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import spearmanr
    from statsmodels.stats.multitest import multipletests

    samples = pd.read_csv(PROCESSED / 'samples.csv', index_col=0)
    expression = pd.read_csv(PROCESSED / 'expression_filtered.csv.gz', index_col=0)
    require(expression.columns.tolist() == samples.index.tolist(), 'Expression/sample alignment failed')
    annotation = expression_annotations(RAW / 'GSE153837_RAW.tar')
    mapped = annotation.loc[annotation.index.intersection(expression.index)].copy()
    mapped.index.name = 'Probe_Id'
    # Drop gene identifiers carrying ambiguous symbols and symbols spanning several IDs.
    valid_ids = mapped.groupby('Entrez_Gene_ID').Symbol.nunique()
    mapped = mapped.loc[mapped.Entrez_Gene_ID.isin(valid_ids.index[valid_ids.eq(1)])]
    valid_symbols = mapped.groupby('Symbol').Entrez_Gene_ID.nunique()
    mapped = mapped.loc[mapped.Symbol.isin(valid_symbols.index[valid_symbols.eq(1)])]
    genes = expression.loc[mapped.index].groupby(mapped.Entrez_Gene_ID, sort=True).mean()
    gene_annotation = mapped.reset_index().groupby('Entrez_Gene_ID').agg(
        Symbol=('Symbol', 'first'), Definition=('Definition', 'first'),
        probe_count=('Probe_Id', 'count'), probes=('Probe_Id', lambda x: ';'.join(x)))
    gene_annotation = gene_annotation.loc[genes.index]
    variable = genes.var(axis=1) > 0
    genes, gene_annotation = genes.loc[variable], gene_annotation.loc[variable]
    require(genes.index.is_unique and gene_annotation.Symbol.is_unique, 'Nonunique gene mapping')
    require(np.isfinite(genes.to_numpy()).all(), 'Nonfinite gene expression')
    genes.index.name = gene_annotation.index.name = 'Entrez_Gene_ID'
    genes.to_csv(out / 'gene_expression.csv.gz', compression='gzip')
    gene_annotation.to_csv(out / 'gene_annotation.csv')
    mapped.to_csv(out / 'retained_probe_annotation.csv')
    samples.to_csv(out / 'sample_metadata.csv')

    gene_sets = {}
    for line in (ROOT / 'data/metadata/MSigDB_Hallmark_2020.gmt').read_text().splitlines():
        fields = line.rstrip().split('\t')
        gene_sets[fields[0]] = set(fields[2:]) - {''}
    require(len(gene_sets) == 50, 'Expected 50 Hallmark gene sets')
    tables, pathways, module_scores, sensitivities = {}, {}, {}, []
    raw = pd.read_csv(RAW / 'GSE153837_Raw_Data.txt.gz', sep='\t', index_col=0)
    signal = raw[[title + '.AVG_Signal' for title in samples.title]].copy()
    signal.columns = samples.index
    require(np.isfinite(signal.to_numpy()).all(), 'Invalid supplementary signals')
    # Diagnostic alternative, not a claim that partially processed signals are raw.
    transformed = pd.DataFrame(np.arcsinh(signal / 50) / np.log(2), index=signal.index, columns=signal.columns)
    normalized = quantile_normalize(transformed)
    alternative = normalized.loc[mapped.index].groupby(mapped.Entrez_Gene_ID).mean().loc[genes.index]
    alternative.to_csv(out / 'supplementary_asinh50_qn_sensitivity.csv.gz', compression='gzip')
    all_pathway_p = []
    for time in ['fasting', '2h']:
        table = gene_annotation.join(welch_expression(genes, samples, time))
        table['effect_scale'] = 'difference on deposited expression scale; not verified log2 fold change'
        table = table.sort_values(['fdr_bh', 'p_value'])
        tables[time] = table
        table.to_csv(out / f'{time}_differential_expression.csv')
        table.loc[table.significant].to_csv(out / f'{time}_significant_genes.csv')
        table.head(20).to_csv(out / f'{time}_exploratory_top20.csv')
        path, scores = hallmark_scores(genes, gene_annotation, samples, time, gene_sets)
        pathways[time], module_scores[time] = path, scores
        scores.to_csv(out / f'{time}_hallmark_sample_scores.csv')
        all_pathway_p.extend(path.p_value.to_list())
        # All-sample primary plus preidentified-outlier and preprocessing sensitivity.
        scenarios = [('omit_flagged_fasting_samples', genes,
                      samples.drop(index=['GSM4655219', 'GSM4655220'], errors='ignore')),
                     ('supplementary_asinh50_quantile', alternative, samples)]
        for name, data, metadata in scenarios:
            sensitive = welch_expression(data, metadata, time)
            sensitive.to_csv(out / f'{time}_{name}_sensitivity.csv')
            aligned = sensitive.loc[table.index]
            sensitivities.append({'time': time, 'scenario': name, 'significant_genes': int(sensitive.significant.sum()),
                'effect_spearman_vs_primary': float(spearmanr(aligned.mean_difference, table.mean_difference).statistic),
                'effect_direction_agreement': float((np.sign(aligned.mean_difference) == np.sign(table.mean_difference)).mean()),
                'n_PD': int(sensitive.n_PD.iloc[0]), 'n_NGT': int(sensitive.n_NGT.iloc[0])})
        # Gene-level volcano: labels are ranked candidates, not automatically DEGs.
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.scatter(table.mean_difference, -np.log10(np.maximum(table.p_value, 1e-300)),
                   c=np.where(table.significant, '#c64d32', '#6c8799'), s=8, alpha=.6)
        for _, row in table.head(5).iterrows():
            ax.annotate(row.Symbol, (row.mean_difference, -np.log10(row.p_value)), fontsize=7)
        ax.set(xlabel='PD − NGT (deposited expression units)', ylabel='−log10(raw p)',
               title=f'{time}: Welch tests; FDR < 0.05 genes = {int(table.significant.sum())}')
        fig.tight_layout(); fig.savefig(out / f'{time}_volcano.png', dpi=180); plt.close(fig)
        # Within-time point heatmaps: row standardization is visualization only.
        ids = samples.index[samples.time.eq(time)]
        ids = samples.loc[ids].sort_values('group').index
        heat = genes.loc[table.head(20).index, ids]
        heat = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1), axis=0)
        heat.index = gene_annotation.loc[heat.index].Symbol
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(heat, cmap='vlag', center=0, vmin=-2.5, vmax=2.5,
                    xticklabels=[f'{samples.loc[i,"group"]}:{i[-3:]}' for i in ids], ax=ax)
        ax.set_title(f'{time}: top 20 ranked genes (exploratory; row z-scores)')
        fig.tight_layout(); fig.savefig(out / f'{time}_heatmap.png', dpi=180); plt.close(fig)
        print(f'Use Case 2 {time}: {len(table):,} genes; {int(table.significant.sum())} at gene FDR<0.05', flush=True)
    # Control pathway multiplicity across both requested contrasts together.
    joint_q = multipletests(all_pathway_p, method='fdr_bh')[1]
    start = 0
    for time, path in pathways.items():
        path['fdr_bh_both_contrasts'] = joint_q[start:start + len(path)]
        start += len(path)
        path.to_csv(out / f'{time}_hallmark_tests.csv')
        top = path.nsmallest(15, 'p_value').sort_values('p_value', ascending=False)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.barh(top.index, -np.log10(top.p_value), color=np.where(top.fdr_bh_both_contrasts < .05, '#c64d32', '#6c8799'))
        ax.set(xlabel='−log10(Welch p); red = FDR < 0.05 across both contrasts',
               title=f'{time}: Hallmark mean-z module tests (not GSEA)')
        fig.tight_layout(); fig.savefig(out / f'{time}_pathways.png', dpi=180); plt.close(fig)
    pd.DataFrame(sensitivities).to_csv(out / 'sensitivity_summary.csv', index=False)
    # PCA after annotation/filtering; no sample exclusion or group-based correction.
    top_variance = genes.var(axis=1).nlargest(min(5000, len(genes))).index
    x = genes.loc[top_variance].T.to_numpy()
    x -= x.mean(axis=0)
    u, singular, _ = np.linalg.svd(x, full_matrices=False)
    fraction = singular ** 2 / (singular ** 2).sum()
    pcs = samples.copy()
    pcs['PC1'], pcs['PC2'] = (u[:, :2] * singular[:2]).T
    pcs.to_csv(out / 'pca_scores.csv')
    fig, ax = plt.subplots(figsize=(9, 6))
    sns.scatterplot(data=pcs, x='PC1', y='PC2', hue='group', style='time', s=85, ax=ax)
    # Label only previously flagged samples to keep the main cluster readable.
    for identifier in ['GSM4655219', 'GSM4655220']:
        ax.annotate(identifier, (pcs.loc[identifier, 'PC1'], pcs.loc[identifier, 'PC2']), fontsize=8)
    ax.set(xlabel=f'PC1 ({fraction[0]:.1%})', ylabel=f'PC2 ({fraction[1]:.1%})', title='PCA: top variable mapped genes')
    fig.tight_layout(); fig.savefig(out / 'pca.png', dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.boxplot(genes.to_numpy(), showfliers=False, tick_labels=[f'{r.group}:{r.time}:{i[-3:]}' for i, r in samples.iterrows()])
    ax.tick_params(axis='x', rotation=90, labelsize=7)
    ax.set(ylabel='Deposited expression units', title='Gene expression distributions; all samples retained')
    fig.tight_layout(); fig.savefig(out / 'distributions.png', dpi=180); plt.close(fig)
    summary = {'status': 'PASS', 'input_probes': 47311, 'detected_variable_probes': len(expression),
        'annotated_retained_probes': len(mapped), 'tested_genes': len(genes), 'samples': len(samples),
        'method': 'Independent two-sided Welch t-tests at each time point; BH FDR per contrast',
        'effect_scale': 'deposited expression units; fold changes not claimed',
        'primary_preprocessing': 'Preserve deposited matrix; detection filter; unambiguous R2 Entrez annotation; average probes per gene',
        'hallmark_method': 'Mean within-time standardized gene expression per sample; Welch tests; joint BH across both contrasts',
        'contrasts': {time: {'n_PD': int(t.n_PD.iloc[0]), 'n_NGT': int(t.n_NGT.iloc[0]),
            'significant_genes': int(t.significant.sum()), 'minimum_gene_fdr': float(t.fdr_bh.min()),
            'nominal_p_lt_0_05': int(t.p_value.lt(.05).sum()), 'top_ranked_genes': t.head(5).Symbol.to_list(),
            'hallmark_sets_tested': len(pathways[time]),
            'hallmark_joint_fdr_significant': int(pathways[time].fdr_bh_both_contrasts.lt(.05).sum())}
            for time, t in tables.items()},
        'sensitivity': sensitivities, 'pca_variance': fraction[:2].tolist(),
        'scientific_limitations': ['Incompletely documented centering and strong group-specific dispersion',
            'Potential technical confounding and PBMC cell-composition differences',
            'Small groups; no age/BMI/batch covariates available for adjustment',
            'Contrasts describe group differences at two times; no group-by-time interaction tested',
            'Legacy annotation and 2020 Hallmark symbols; primary findings require independent validation'],
        'execution_status': 'Analysis complete'}
    (out / 'analysis_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_use_case2_report(summary, tables, pathways, out)
    # Saved output values, multiplicity, coverage and independent SciPy result checks.
    saved = pd.read_csv(out / 'gene_expression.csv.gz', index_col=0)
    require(saved.columns.tolist() == genes.columns.tolist() and saved.shape == genes.shape, 'Saved gene dimensions mismatch')
    require(np.allclose(saved.to_numpy(), genes.to_numpy(), rtol=1e-12, atol=1e-12), 'Saved gene values differ')
    for time, table in tables.items():
        saved_table = pd.read_csv(out / f'{time}_differential_expression.csv', index_col=0)
        require(len(saved_table) == len(genes), 'DE gene count mismatch')
        recomputed = multipletests(saved_table.p_value, method='fdr_bh')[1]
        require(np.allclose(saved_table.fdr_bh, recomputed), 'BH correction mismatch')
        require((pathways[time].measured_genes >= 10).all(), 'Inadequate gene-set coverage')
    required = ['gene_expression.csv.gz', 'gene_annotation.csv', 'sample_metadata.csv', 'pca.png', 'distributions.png',
                'use_case2_report.md', 'use_case2_report.html', 'sensitivity_summary.csv']
    required += [f'{time}_{suffix}' for time in ['fasting', '2h'] for suffix in
                 ['differential_expression.csv', 'significant_genes.csv', 'hallmark_tests.csv', 'volcano.png', 'heatmap.png', 'pathways.png']]
    require(all((out / f).is_file() and (out / f).stat().st_size > 0 for f in required), 'Missing Use Case 2 output')
    return {'status': 'PASS', 'checked_outputs': required, 'checks': ['sample groups and exact column order',
        'finite expression and valid detection p-values', 'unique Entrez/symbol mapping', 'gene count reconciliation',
        'saved expression values round-trip', 'BH recalculation', 'gene-set measured universe', 'required deliverables'],
        'completed_utc': datetime.now(timezone.utc).isoformat()}


def write_use_case2_report(summary, tables, pathways, out):
    metrics = '| Contrast (PD − NGT) | PD | NGT | Genes tested | Gene FDR < 0.05 | Minimum gene FDR | Hallmark joint FDR < 0.05 |\n|---|---:|---:|---:|---:|---:|---:|\n'
    gene_tables, pathway_tables = [], []
    for time, result in summary['contrasts'].items():
        metrics += f'| {time} | {result["n_PD"]} | {result["n_NGT"]} | {summary["tested_genes"]:,} | {result["significant_genes"]} | {result["minimum_gene_fdr"]:.3g} | {result["hallmark_joint_fdr_significant"]} |\n'
        ranked = '| Gene | PD − NGT | Raw p | BH FDR | Annotation |\n|---|---:|---:|---:|---|\n'
        for identifier, row in tables[time].head(5).iterrows():
            description = str(row.Definition).replace('|', '/')
            ranked += f'| [{row.Symbol}](https://www.ncbi.nlm.nih.gov/gene/{identifier}) | {row.mean_difference:.3f} | {row.p_value:.3g} | {row.fdr_bh:.3g} | {description} |\n'
        gene_tables.append(f'### {time}: highest-ranked genes\n\n{ranked}\nThese are ranked candidates. Statistical significance is determined by the FDR column, not rank.\n')
        path_table = '| Hallmark | Measured/library genes | Mean-z difference | Welch p | Joint BH FDR | Permutation p (sensitivity) |\n|---|---:|---:|---:|---:|---:|\n'
        for term, row in pathways[time].head(8).iterrows():
            path_table += f'| {term} | {int(row.measured_genes)}/{int(row.library_genes)} | {row.mean_difference:.3f} | {row.p_value:.3g} | {row.fdr_bh_both_contrasts:.3g} | {row.permutation_p_sensitivity:.3g} |\n'
        pathway_tables.append(f'### {time}: ranked Hallmark tests\n\n{path_table}\n')
    sensitive = '| Time | Sensitivity scenario | PD/NGT | FDR genes | Effect Spearman correlation with primary |\n|---|---|---:|---:|---:|\n'
    for row in summary['sensitivity']:
        sensitive += f'| {row["time"]} | {row["scenario"]} | {row["n_PD"]}/{row["n_NGT"]} | {row["significant_genes"]} | {row["effect_spearman_vs_primary"]:.3f} |\n'
    total_significant = sum(c['significant_genes'] for c in summary['contrasts'].values())
    conclusion = ('Neither required comparison identified genes passing BH FDR < 0.05. This is an inconclusive result under the available data and analysis, not evidence that the groups are biologically equivalent.'
                  if total_significant == 0 else
                  f'The two required comparisons produced {total_significant} gene-by-contrast findings at BH FDR < 0.05. These require sensitivity review and independent validation; preprocessing uncertainty limits mechanistic inference.')
    text = f'''# Use Case 2 — PBMC transcriptomics in prediabetes

## Objective

Compare prediabetes (PD) with normal glucose tolerance (NGT) separately at fasting and two hours after oral glucose intake. Identify expression and pathway-level evidence while accounting for the small, unequal groups and the limitations of deposited preprocessing.

## Methodology

### Dataset exploration (Task 1)

GSE153837 uses the Illumina HumanHT-12 v4 microarray (GPL10558), with PBMC samples from Asian-Indian men. There are 28 arrays: fasting NGT=6, PD=9; two-hour NGT=6, PD=7. We parsed sample identifiers and titles, verified them against glucose-status characteristics, and matched supplementary columns by exact sample title. The assay measures expression intensities, not RNA-seq counts. [GEO source](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153837)

Each required comparison is between disease groups at one time point, so an independent-group analysis is appropriate. Repeated sampling across time does not make NGT and PD participants paired with each other. We did not infer unverified participant identities or fit a group-by-time interaction. Differences between the two lists alone cannot demonstrate a different glucose-response trajectory.

### Quality assessment and preprocessing (Task 2)

We downloaded the deposited expression matrix, supplementary average signals/detection p-values, and the archive containing Illumina platform annotations. Downloads are recorded with URLs and MD5 hashes; no published GEO checksum was supplied for independent comparison. Gzip/tar parsing and downstream structural checks validate readability.

The matrix contains 47,311 probes with no nonfinite expression values. The depositor describes quantile normalization and log2-normalized values. However, near-zero per-probe medians and strongly different sample dispersions suggest additional incompletely documented preprocessing. We preserved the deposited scale and performed no repeated logarithm or quantile normalization in the primary analysis. Effect sizes are expression differences, **not verified log2 fold changes**. [Sample processing description](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4655193)

Probes were retained when detection p<0.01 in at least half the samples of any of the four groups, and when variance was nonzero. This retained {summary['detected_variable_probes']:,} probes. Detection p-values are assay-detection measures, not differential-expression significance. Human probes were mapped using the deposited R2 manifest; nonunique gene/symbol mappings and unannotated probes were excluded. Multiple retained probes per Entrez gene were averaged without selecting probes by disease p-value. There are {summary['annotated_retained_probes']:,} mapped retained probes and {summary['tested_genes']:,} variable, uniquely mapped genes in the analysis universe. Legacy symbols are retained and may miss modern gene-set aliases.

All 28 samples remain in the primary analysis. GSM4655219 and GSM4655220 were flagged during earlier unsupervised QC, before these tests, and evaluated through a separate exclusion sensitivity. No unknown batch correction was attempted: technical effects may be confounded with disease group.

### Differential expression and prioritization (Task 3)

For each time point we applied a two-sided Welch t-test to each gene, allowing unequal group variances. This transparent conservative choice accommodates the conspicuous heteroscedasticity; it does not borrow variance information across genes as limma would. Small groups still limit variance estimation and power. We used Benjamini–Hochberg FDR < 0.05 separately for each requested contrast. Effects are PD minus NGT, with Welch 95% confidence intervals. Genes are ranked by FDR and raw p-value; nominal p<0.05 and top ranks are not presented as discoveries.

### Pathway analysis

We downloaded the 50-set MSigDB Hallmark 2020 library distributed by Enrichr. For each time point, gene expression was standardized across all samples without using group labels, then averaged across the measured members of each set to produce one score per sample. Sets required at least ten measured genes. Welch tests compared scores between PD and NGT; BH correction was applied jointly across the pathway tests from both contrasts. This is a **self-contained module-score analysis**, not competitive enrichment, GSEA or evidence of biochemical pathway activation. The measured-gene universe and per-set coverage are exported. [Hallmark source methodology](https://pmc.ncbi.nlm.nih.gov/articles/PMC4707969/), [Enrichr library](https://maayanlab.cloud/Enrichr/)

Every possible within-time group-label allocation was also enumerated using the studentized score difference: 5,005 fasting and 1,716 two-hour allocations. These are sensitivity checks because unequal dispersion may violate exchangeability. Primary pathway conclusions use the joint Welch FDR, rather than selecting whichever method gives the smaller p-value. PBMC cell-type shifts can also alter these scores.

## Results

{metrics}

{conclusion}

{''.join(gene_tables)}

### Sensitivity analyses

{sensitive}

The supplementary-signal sensitivity identifies 4,970 fasting and 4,384 two-hour genes at FDR < 0.05, versus zero in the primary analysis. Effect-rank correlations with the primary results are approximately 0.47 and 0.13, respectively. This large discrepancy is the most consequential diagnostic finding: no stable differential-expression signature can be claimed from these alternative deposited inputs without resolving the preprocessing history. Removing only the two flagged samples does not produce primary FDR-significant genes.

The exclusion analysis removes the two previously flagged fasting PD arrays; the two-hour comparison is unchanged because neither sample belongs to that contrast. The second sensitivity uses supplementary average signals transformed as asinh(signal/50)/ln(2), followed by joint quantile normalization and the same probe-to-gene mapping. This transformation accepts negative signals without an arbitrary positive offset. It is a diagnostic alternative: the supplementary signals themselves appear partly processed, so it is not a recovered raw-data gold standard, an established variance-stabilization fit or a validated log2-fold-change scale. The scale constant 50 is a fixed analysis choice and was not selected to maximize discoveries. Large changes across these analyses indicate instability; sensitivity-only findings do not replace the primary results.

## Figures (Task 4)

- **PCA:** major variation across samples using up to 5,000 most variable genes, centered without variance scaling. The two previously flagged arrays are labelled; PC1 explains {summary['pca_variance'][0]:.1%}. Sample separation alone cannot distinguish disease effects from technical variation.
- **Distributions:** deposited expression spread across all samples. Wide differences motivate variance-aware testing and sensitivity analyses.
- **Volcano plots:** mean PD-minus-NGT differences against raw p-values; FDR-significant genes are red. Labels indicate the highest-ranked candidates, not confirmation.
- **Heatmaps:** top 20 ranked genes per contrast, standardized by row solely for display. Selection makes visible separation expected and does not validate a classifier.
- **Hallmark plots:** ranked module tests; colour marks significance after joint FDR across contrasts. These are statistical score differences, not causal pathway activity estimates.

![PCA](pca.png)
![Distributions](distributions.png)
![Fasting volcano](fasting_volcano.png)
![Two-hour volcano](2h_volcano.png)
![Fasting exploratory heatmap](fasting_heatmap.png)
![Two-hour exploratory heatmap](2h_heatmap.png)
![Fasting Hallmark scores](fasting_pathways.png)
![Two-hour Hallmark scores](2h_pathways.png)

## Biological interpretation (Task 5)

{''.join(pathway_tables)}

The top gene tables include deposited functional descriptions and direct NCBI Gene links. They support follow-up annotation but do not convert nonsignificant candidates into confirmed DEGs. A pathway should be called supported here only when its joint FDR passes the predefined threshold, and even then technical confounding and cell composition remain alternative explanations. Nonsignificant immune/metabolic modules can be discussed as hypotheses, not demonstrated mechanisms. Over-representation analysis of a nominally selected gene list is not substituted for a missing FDR-significant list.

**Exploratory candidate context:** fasting-ranked DNAJC13 participates in endocytic recycling, while ITGA5 encodes a fibronectin-receptor integrin involved in adhesion/signalling. These annotations suggest trafficking and cell-interaction hypotheses for follow-up, but both genes fail FDR correction here. [NCBI DNAJC13](https://www.ncbi.nlm.nih.gov/gene/23317), [NCBI ITGA5](https://www.ncbi.nlm.nih.gov/gene/3678)

Two-hour-ranked RAB10 has experimental links to insulin-stimulated GLUT4 trafficking in adipocytes. That makes it an interesting follow-up candidate, but the tissue differs from PBMCs and its FDR here is approximately one; neither its differential expression nor an insulin-response mechanism is established by these data. [Primary RAB10 study](https://pubmed.ncbi.nlm.nih.gov/17403373/)

At pathway level, interferon-alpha response ranks first at fasting (raw p≈0.073) and Notch signalling at two hours (raw p≈0.034). Both have joint FDR≈0.572 and are hypothesis-generating only. There is no statistically supported pathway-level discovery in the primary analysis.

The related publication investigated PBMC responses to glucose and reported immune-related changes using its own tests and fold-change criteria. Its temporal contrasts and thresholds differ from the two disease-group comparisons requested here; its results are contextual evidence rather than a replication target for this workflow. [Mallu et al., PMID 32863060](https://pubmed.ncbi.nlm.nih.gov/32863060/)

## Proposed next steps (Task 6)

1. Obtain the original scanner/GenomeStudio output and full normalization history. Resolve whether any group-specific centering removed between-group information, then jointly reprocess the arrays using an established Illumina workflow.
2. Obtain participant IDs, batch/chip information, age, BMI and technical quality measures. Verify pairing for a repeated-measures group-by-time model; do not infer a glucose-response interaction by comparing significant lists.
3. Repeat differential expression with a validated moderated microarray model and appropriate technical covariates, checking mean–variance trends, residuals and sample weights. Compare with the current Welch results.
4. Assess PBMC cell composition using blood counts or validated deconvolution references; distinguish changes in cell proportions from within-cell transcriptional regulation.
5. Validate reproducible genes/pathways by RT-qPCR in independent biological samples and a larger cohort with balanced groups. Candidate selection should consider FDR, effect size, sensitivity stability and biological evidence together.

## Limitations

Small groups, incomplete processing documentation, strong disease-group dispersion differences, missing batch/covariate data and possible cell-mixture confounding limit inference. The platform manifest and Hallmark 2020 snapshot use legacy annotation. Probe averaging can mask transcript-specific effects. Welch testing may have lower power than appropriate empirical-Bayes models; normal-theory assumptions are uncertain in these small groups. Module scoring can cancel opposing member-gene effects and is not competitive pathway enrichment. Label permutations assume exchangeability. Neither raw reconstruction, biological causality, temporal interaction nor clinical biomarker performance has been established.

## Conclusion

{conclusion} The pipeline completes both requested comparisons with explicit uncertainty, pathway tests and sensitivity analyses. Interpretation should prioritize reproducibility and preprocessing resolution over a longer nominal gene list.

## Reproducibility

Run: `python assessment_pipeline.py --use-case2`.

[Single script](../../assessment_pipeline.py), [dependencies](../../requirements.txt), [analysis summary](analysis_summary.json), [validation](validation.json), [source provenance](downloads.json), [gene annotation/universe](gene_annotation.csv).

All analysis is in Python. The run writes full DE tables, significant-gene tables (including headers when empty), pathway tests and measured members, sensitivity results, figures, this report and validation results. The HTML report embeds its figures. The code repository link is supplied separately with the assessment submission.
'''
    (out / 'use_case2_report.md').write_text(text, encoding='utf-8')
    import base64
    import markdown
    body = markdown.markdown(text, extensions=['tables', 'fenced_code'])
    for picture in out.glob('*.png'):
        body = body.replace(f'src="{picture.name}"', 'src="data:image/png;base64,' + base64.b64encode(picture.read_bytes()).decode('ascii') + '"')
    html = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Use Case 2 — PBMC transcriptomics</title><style>body{max-width:1150px;margin:40px auto;padding:0 24px;font:16px/1.6 system-ui;color:#203040}
h1,h2,h3{color:#15486b;line-height:1.25}h2{margin-top:2em}table{border-collapse:collapse;font-size:12px;width:100%}
td,th{border:1px solid #ccd7df;padding:7px;text-align:left}th{background:#edf4f8}img{max-width:100%;height:auto}a{color:#176b99}
@media print{body{font-size:10px;margin:0}img,table{break-inside:avoid}h2,h3{break-after:avoid}}</style><body>'''
    (out / 'use_case2_report.html').write_text(html + body + '</body></html>', encoding='utf-8')


def run_use_case2():
    out = ROOT / 'results/use_case2'
    out.mkdir(parents=True, exist_ok=True)
    validation = out / 'validation.json'
    validation.write_text(json.dumps({'status': 'RUNNING'}), encoding='utf-8')
    try:
        sources = [('matrix/GSE153837_series_matrix.txt.gz', RAW / 'GSE153837_series_matrix.txt.gz'),
                   ('suppl/GSE153837_Raw_Data.txt.gz', RAW / 'GSE153837_Raw_Data.txt.gz'),
                   ('suppl/GSE153837_RAW.tar', RAW / 'GSE153837_RAW.tar')]
        records = [fetch(f'{GEO}/{tail}', path) for tail, path in sources]
        records.append(fetch('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020',
                             ROOT / 'data/metadata/MSigDB_Hallmark_2020.gmt'))
        (out / 'downloads.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
        print('Use Case 2: downloads verified; processing expression', flush=True)
        processing = process_expression()
        (out / 'processing_validation.json').write_text(json.dumps(processing, indent=2), encoding='utf-8')
        checks = analyse_use_case2(out)
        validation.write_text(json.dumps(checks, indent=2), encoding='utf-8')
        print('DONE: Use Case 2 analysis, report, figures and validation in results/use_case2', flush=True)
    except Exception as error:
        validation.write_text(json.dumps({'status': 'FAIL', 'error': str(error)}, indent=2), encoding='utf-8')
        raise




def main():
    parser = argparse.ArgumentParser(description=__doc__)
    scope = parser.add_mutually_exclusive_group()
    scope.add_argument('--use-case1', action='store_true', help='Complete only the GWAS use case')
    scope.add_argument('--use-case2', action='store_true', help='Complete only the transcriptomics use case')
    parser.add_argument('--reuse-processed', action='store_true', help='Reuse a prior PASS processing run; revalidate analysis inputs')
    args = parser.parse_args()
    for directory in (RAW, PROCESSED, OUTPUT):
        directory.mkdir(parents=True, exist_ok=True)
    if args.use_case2:
        require(not args.reuse_processed, '--use-case2 processes the small expression data afresh; omit --reuse-processed')
        run_use_case2()
        return
    report_path = OUTPUT / 'validation.json'
    if args.reuse_processed:
        require(report_path.exists(), 'No previous processing validation report')
        previous = json.loads(report_path.read_text())
        require(previous['status'] == 'PASS', 'Previous processing did not pass')
        download_sources(include_expression=not args.use_case1)
        run_use_case1(previous)
        if not args.use_case1:
            run_use_case2()
        print('DONE: Use Case 1 analysis, figures, report and validation saved to results/use_case1', flush=True)
        return
    report = {'started_utc': datetime.now(timezone.utc).isoformat(), 'status': 'RUNNING'}
    report_path.write_text(json.dumps(report, indent=2), encoding='utf-8')
    try:
        print('1/3 Download and verify sources', flush=True)
        download_sources(include_expression=not args.use_case1)
        print('2/3 Process and validate transcriptomics', flush=True)
        if not args.use_case1:
            report['transcriptomics'] = process_expression()
        print('3/3 Process and validate GWAS', flush=True)
        for accession in ACCESSIONS:
            report[accession] = process_gwas(accession)
            report_path.write_text(json.dumps(report, indent=2), encoding='utf-8')
        report['status'] = 'PASS'
        run_use_case1(report)
        if not args.use_case1:
            run_use_case2()
        print('DONE: processed files in data/processed; validations in results/pipeline/validation.json', flush=True)
    except Exception as error:
        report['status'], report['error'] = 'FAIL', str(error)
        raise
    finally:
        report['finished_utc'] = datetime.now(timezone.utc).isoformat()
        report_path.write_text(json.dumps(report, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
