# CRC Single-Cell Atlas: Full Condition-Level Integration

This repository contains the condition-level integration workflow used for the full colorectal cancer (CRC) single-cell atlas.

The workflow integrates the post-QC study-condition Seurat objects separately within four biological conditions:

- `Adenoma_polyp`
- `CRC`
- `Normal`
- `TAN` (tumor-adjacent normal)

`Serrated_polyp` is not cross-study integrated because only one study contributes to that condition.

## Scope

This workflow is for the **full post-QC atlas**, not the 20,000-cell subsampled integration.

Expected full-atlas cell counts are:

| Condition | Contributing studies | Post-QC cells |
|---|---:|---:|
| Adenoma_polyp | 3 | 306,396 |
| CRC | 13 | 702,029 |
| Normal | 5 | 151,684 |
| TAN | 7 | 287,894 |
| Serrated_polyp | 1 | 18,779 |

The four multi-study conditions are integrated with both:

- Seurat CCA anchor-based integration
- Harmony PCA-space batch correction

## Input data

The workflow expects post-QC study-condition Seurat objects under:

```text
<atlas_root>/full_seurat_QC_normalized_by_study_condition/
```

with filenames:

```text
<Study>_<Condition>_QC_normalized_seurat.rds
```

Examples:

```text
GSE201348_Adenoma_polyp_QC_normalized_seurat.rds
GSE132465_CRC_QC_normalized_seurat.rds
GSE294300_TAN_QC_normalized_seurat.rds
```

Condition-specific HVG files are expected under:

```text
<atlas_root>/full_atlas_HVG_2500/
```

with filenames:

```text
Adenoma_polyp_2500_HVGs.csv
CRC_2500_HVGs.csv
Normal_2500_HVGs.csv
TAN_2500_HVGs.csv
```

Each HVG file must contain a column called:

```text
gene
```

Up to 2,500 condition-specific HVGs are supplied. Only HVGs present in the expression representation of **every contributing study within that condition** are retained for integration.

## Study composition

### Adenoma_polyp

- `CELLxGENE_Chen_2021`
- `GSE201348`
- `GSE261388`

### CRC

- `CELLxGENE_Chen_2021`
- `CELLxGENE_Moorman2025`
- `GSE132465`
- `GSE178341`
- `GSE200997`
- `GSE201348`
- `GSE221575`
- `GSE225857`
- `GSE231559`
- `GSE261388`
- `GSE294300`
- `GSE302903`
- `GSE330797`

### Normal

- `CELLxGENE_Chen_2021`
- `GSE201348`
- `GSE221575`
- `GSE231559`
- `GSE261388`

### TAN

- `CELLxGENE_Moorman2025`
- `GSE132465`
- `GSE178341`
- `GSE200997`
- `GSE225857`
- `GSE231559`
- `GSE294300`

## Preprocessing rules preserved during integration

The input cohorts were not all supplied on the same preprocessing scale, so the integration workflow preserves their original provenance rather than renormalizing everything indiscriminately.

### GEO-derived datasets

The GEO-derived study-condition objects were previously log-normalized.

The workflow therefore uses the existing normalized RNA `data` layers and does **not** call `NormalizeData()` again.

### CELLxGENE Chen 2021

The Chen cohort was supplied as an already processed/scaled matrix stored in the RNA `counts` layer.

Therefore:

- `NormalizeData()` is not rerun.
- The supplied Chen expression matrix is not rerun through `ScaleData()`.
- For CCA, the supplied processed Chen matrix is used as the study expression representation.
- For Harmony, PCA is learned from the non-Chen datasets and the Chen matrix is projected into that PCA space using the non-Chen PCA loadings.

The projection used is conceptually:

```r
chen_pca <- Matrix::t(chen_scaled) %*% pca_loadings
```

### CELLxGENE Moorman 2025

The Moorman cohort was already normalized.

Therefore:

- `NormalizeData()` is not rerun.
- If a normalized RNA `data` layer is available, it is used.
- If normalized values are stored in the RNA `counts` layer, that layer is used as the normalized expression representation.

This is an integration/visualization choice only. A normalized Moorman matrix should not be treated as raw integer counts for DESeq2.

## Integration features

For each condition:

1. Load the condition-specific HVG list.
2. Determine which HVGs are present in every contributing study.
3. Retain the shared HVGs only.
4. Use the same shared feature set for CCA and Harmony.

The final feature list is written to:

```text
<Condition>_shared_integration_HVGs.csv
```

## CCA workflow

CCA integration uses Seurat's anchor-based workflow:

```text
study-specific normalized/processed expression
        |
        v
shared condition-specific HVGs
        |
        v
FindIntegrationAnchors(reduction = "cca")
        |
        v
IntegrateData()
        |
        v
integrated assay
        |
        v
ScaleData()
        |
        v
PCA (50 PCs)
        |
        v
UMAP / neighbors using dimensions 1:30
        |
        v
clustering at resolutions 0.3 and 0.5
```

The CCA reductions and graph names are:

```text
pca.cca
umap.cca
cca_nn
cca_snn
```

Cluster columns are:

```text
cca_snn_res.0.3
cca_snn_res.0.5
```

The CCA integrated assay contains corrected integration values. These values are **not raw RNA counts** and should not be used as raw-count differential-expression input.

## Harmony workflow

Harmony corrects the PCA representation rather than generating a corrected gene-expression assay.

For conditions without Chen:

```text
normalized RNA expression
        |
        v
ScaleData(shared HVGs)
        |
        v
PCA (50 PCs)
        |
        v
Harmony by study
        |
        v
UMAP / neighbors using Harmony dimensions 1:30
        |
        v
clustering at resolutions 0.3 and 0.5
```

For conditions containing Chen:

```text
non-Chen normalized RNA expression
        |
        v
ScaleData(shared HVGs)
        |
        v
PCA (50 PCs)
        |
        +--------------------------+
        |                          |
        |                          v
        |               Chen processed/scaled matrix
        |                          |
        |                          v
        |               project using PCA loadings
        |                          |
        +-------------+------------+
                      |
                      v
          combined PCA representation
                      |
                      v
              Harmony by study
                      |
                      v
         UMAP / neighbors / clustering
```

Harmony uses:

```text
reduction: harmony
UMAP: umap.harmony
graphs: harmony_nn, harmony_snn
```

Cluster columns are:

```text
harmony_snn_res.0.3
harmony_snn_res.0.5
```

The Harmony object is primarily an integrated dimensional representation. It should not be treated as a corrected gene-by-cell expression matrix for differential expression.

## Dimensionality reduction and clustering

For both integration approaches:

- 50 principal components are calculated.
- Dimensions 1 through 30 are used downstream.
- UMAP is calculated from those dimensions.
- A nearest-neighbor graph is constructed.
- Clustering is evaluated at resolutions 0.3 and 0.5.

## Output structure

The script writes results under:

```text
<atlas_root>/Integration_by_condition/
```

For example:

```text
Integration_by_condition/
├── Adenoma_polyp/
│   ├── Adenoma_polyp_shared_integration_HVGs.csv
│   ├── Adenoma_polyp_CCA_integrated.rds
│   └── Adenoma_polyp_Harmony_integrated.rds
├── CRC/
│   ├── CRC_shared_integration_HVGs.csv
│   ├── CRC_CCA_integrated.rds
│   └── CRC_Harmony_integrated.rds
├── Normal/
│   ├── Normal_shared_integration_HVGs.csv
│   ├── Normal_CCA_integrated.rds
│   └── Normal_Harmony_integrated.rds
└── TAN/
    ├── TAN_shared_integration_HVGs.csv
    ├── TAN_CCA_integrated.rds
    └── TAN_Harmony_integrated.rds
```

The workflow also writes:

```text
full_atlas_integration_manifest.csv
full_atlas_integration_sessionInfo.txt
```

## Differential expression

The integrated objects are intended for:

- batch-corrected visualization
- neighborhood construction
- UMAP
- clustering
- transfer of cluster/cell-state assignments

Differential expression should be performed using the appropriate source RNA expression or raw-count matrices rather than the CCA integrated assay.

For the atlas pseudobulk DESeq2 analysis, CCA cluster assignments were transferred back to the study-specific raw-count objects by exact cell identifiers.

The Chen and Moorman atlas matrices were excluded from raw-count DESeq2 pseudobulk because the versions used in this atlas were already processed/normalized rather than valid raw integer count inputs.

## Requirements

R packages:

```r
Seurat
SeuratObject
Matrix
data.table
harmony
```

The workflow uses Seurat v5-style assay layers.

## Running the workflow

On OSCER:

```bash
export CRC_ATLAS_ROOT=/ourdisk/hpc/moussalab/jumi/dont_archive/atlas
Rscript 01_full_condition_integration.R
```

If `CRC_ATLAS_ROOT` is not set, the script defaults to the atlas path used on OSCER.

The analysis is computationally intensive because the full CRC condition contains more than 700,000 cells. Run it through an appropriate OSCER/Slurm allocation rather than a login node.

## Configuration

At the top of `01_full_condition_integration.R`, the following options can be changed:

```r
config$run_cca
config$run_harmony
config$save_intermediate
config$overwrite
```

Setting:

```r
config$overwrite <- FALSE
```

prevents existing final integration files from being overwritten.

## Reproducibility note

This repository version is a cleaned, reproducible implementation of the documented full-atlas integration workflow. It preserves the atlas-specific preprocessing decisions, study composition, HVG strategy, dimensionality, and clustering parameters used in the project.

Because the historical analysis was developed interactively across multiple condition-specific runs and repairs, this repository script should be treated as the consolidated workflow rather than a byte-for-byte copy of every original interactive command.

## Recommended repository layout

```text
crc-single-cell-atlas-integration/
├── README.md
└── 01_full_condition_integration.R
```

## Data availability

Large Seurat RDS files are not intended to be committed directly to GitHub. The repository should contain the workflow code, documentation, manifests, and small metadata tables, while large input/output objects remain on institutional storage or an appropriate data repository.
