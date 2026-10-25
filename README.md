# mohan-s-assessment
This repository contains the code base for MDRF assessment 

# Bioinformatics assessment

Reproducible Python workflow for GWAS summary statistics and PBMC microarray analysis.
Validated with Python 3.12.3. Internet is required for initial source downloads.
All paths are relative to the script location. No existing data/output folders are required.
Allow several GB of free disk space; the two GWAS downloads total approximately 845 MB.

## Setup

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

On Linux/macOS, use `.venv/bin/python` instead of `.\.venv\Scripts\python.exe`.

## Execute

Both use cases:

```powershell
.\.venv\Scripts\python.exe assessment_pipeline.py
```

Or run either independently:

```powershell
.\.venv\Scripts\python.exe assessment_pipeline.py --use-case1
.\.venv\Scripts\python.exe assessment_pipeline.py --use-case2
```

The script downloads public inputs, processes data, performs analyses, writes reports,
and runs output validations. Failures return a nonzero exit code and are recorded in
the relevant validation JSON. Completed downloads are reused; interrupted transfers resume.
No ZIP archives are generated.

To regenerate GWAS analysis from an existing successful processing run:

```powershell
.\.venv\Scripts\python.exe assessment_pipeline.py --use-case1 --reuse-processed
```

Do not use `--reuse-processed` for a fresh checkout.

## Outputs and validation

| Location | Contents |
|---|---|
| `data/raw/` | Downloaded source files |
| `data/metadata/` | Source metadata/checksums and gene-set definitions |
| `data/processed/` | Cleaned GWAS and filtered expression data |
| `results/pipeline/` | Processing counts, download manifest and GWAS processing validation |
| `results/use_case1/` | GWAS tables, seven figures, scientific report and analysis validation |
| `results/use_case2/` | Expression tables, eight figures, scientific report and analysis validation |

Open `use_case1_report.html` and `use_case2_report.html` in their respective output folders.
Each report embeds its figures. Inspect `validation.json` in each use-case folder:
`status` should be `PASS`. PASS confirms technical checks, not biological significance.

### Reference results

| Check | Expected result for the supplied source snapshots |
|---|---|
| Primary GWAS retained variants | 8,340,560 |
| BMI-adjusted GWAS retained variants | 8,338,623 |
| Primary/BMI-adjusted significant variants (p < 5e-8) | 4 / 7 |
| Descriptive genomic regions | 2, near GCK and HK1 |
| Expression samples / filtered probes / mapped genes | 28 / 17,652 / 12,883 |
| Fasting comparison: PD / NGT | 9 / 6 |
| Two-hour comparison: PD / NGT | 7 / 6 |
| Primary genes at BH FDR < 0.05 | 0 in both contrasts |
| Hallmark modules at joint FDR < 0.05 | 0 in both contrasts |

Compare GWAS counts and inflation statistics in `results/use_case1/analysis_summary.json`;
compare expression counts and sensitivity results in `results/use_case2/analysis_summary.json`.
Floating-point values may differ slightly between platforms. Source manifests record input hashes.

## Methods and interpretation

GWAS: published MD5 verification, chunked QC, p < 5e-8 prioritization, odds ratios/intervals,
allele-aligned BMI sensitivity analysis, Manhattan/Q-Q/regional/effect plots. Distance grouping
is descriptive; LD independence and causality are not inferred. Models share participants.

Transcriptomics: preserve deposited expression scale; detection filtering; unambiguous
R2-manifest Entrez mapping; average probes per gene; independent two-sided Welch tests
with BH correction at each time point. Pathway analyses use sample-level Hallmark mean-z
scores with joint BH correction across contrasts; they are not GSEA or competitive enrichment.
Sensitivity analyses assess two flagged samples and supplementary-signal preprocessing.
The latter changes findings substantially, so no stable expression signature is claimed.
Effects on the deposited scale are not interpreted as verified log2 fold changes.

Biological interpretation in the generated reports is curated and cited, rather than
an automated gene-function discovery step. No genotype-level GWAS QC, temporal interaction
model, independent replication or clinical prediction is claimed.

## Submission

The code repository needs only `assessment_pipeline.py`, `requirements.txt`, `README.md`
and `.gitignore`. Deliver the two generated HTML scientific reports separately as assessment
reports. Downloaded data, generated outputs, environments and exploratory scripts remain local.
