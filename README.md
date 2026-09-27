cat > README.md <<'EOF'
# Wilms Tumor RNA-seq Analysis

R scripts for differential gene expression and pathway enrichment analysis of Wilms tumor RNA-seq data.

The analysis compares:
- Primary tumors (PRIM) vs. normal tissue (NORM)
- Recurrent tumors (RECU) vs. normal tissue (NORM)

The repository preserves the original research analysis and provides instructions for reusing the code with a compatible gene-count matrix and sample metadata.

## Repository Contents

| File | Description |
|---|---|
| `merged_counts.csv` | Gene-count matrix used as input |
| `sample_info.tsv` | Sample metadata |
| `install.R` | Package installation script |
| `deseq2_script.R` | Primary tumor vs. normal differential expression analysis |
| `deseq2_script_recu.R` | Recurrent tumor vs. normal differential expression analysis |
| `pathway_all.R` | Supplementary pathway analysis across comparisons |
| `pathway_gsea.R` | Supplementary pathway enrichment and visualization code |

## Requirements

- R
- The R packages specified in `install.R`

Install dependencies from the repository's root directory:

```bash
Rscript install.R
```

## Input Data

The analysis expects two files in the repository's root directory.

### Gene counts: `merged_counts.csv`

A CSV file containing raw gene counts:

- Rows represent genes.
- The first column contains gene identifiers.
- Remaining columns contain sample counts.
- Sample column names must match the sample identifiers in `sample_info.tsv`.

The scripts use Ensembl gene identifiers for downstream annotation and enrichment analysis.

### Sample metadata: `sample_info.tsv`

A tab-separated file containing sample information.

The original scripts expect a `sample` column and a `condition` column. The condition values used by the code are:

- `NORM`: normal tissue
- `PRIM`: primary tumor
- `RECU`: recurrent tumor

If using a different dataset, keep the input structure and condition labels compatible with the scripts, or update the corresponding references in the code.

## Running the Analysis

Clone the repository and enter its directory:

```bash
git clone <repository-url>
cd wilms-tumor-analysis
```

Install the dependencies:

```bash
Rscript install.R
```

Run the primary tumor vs. normal analysis:

```bash
Rscript deseq2_script.R > primary_run.log 2>&1
```

Run the recurrent tumor vs. normal analysis:

```bash
Rscript deseq2_script_recu.R > recurrent_run.log 2>&1
```

These commands should be run from the repository's root directory because the scripts use relative paths to locate their input files.

The log files capture analysis progress and verbose output. They are excluded from Git by `.gitignore`.

## Analysis Overview

The DESeq2 scripts:

1. Load the raw count matrix and sample metadata.
2. Select samples for the relevant comparison.
3. Construct a DESeq2 dataset using `condition` as the design variable.
4. Filter genes with low counts.
5. Estimate size factors and run differential expression analysis.
6. Generate exploratory visualizations and downstream analysis outputs.

The scripts use `NORM` as the reference condition.

## Supplementary Pathway Scripts

`pathway_all.R` and `pathway_gsea.R` contain additional pathway enrichment and visualization code from the research workflow.

These scripts currently expect DESeq2 result objects to exist in the R environment. They are included as supplementary research code and are not standalone entry points for a fresh R session.

## Reusing the Code

To analyze another compatible dataset:

1. Replace `merged_counts.csv` with your raw count matrix.
2. Replace `sample_info.tsv` with matching sample metadata.
3. Preserve the expected sample identifiers and condition labels, or modify the scripts accordingly.
4. Run the relevant DESeq2 script from the repository's root directory.

Do not use normalized expression values in place of raw counts as input to `DESeqDataSetFromMatrix()`.

## Notes

This repository contains research analysis scripts rather than a packaged R library or automated production pipeline. Some supplementary visualizations were originally developed interactively in RStudio.
EOF