# SpaRGS: Spatial Rank-based Gene-set Scoring

This repository contains the SpaRGS implementation and analysis workflows supporting the manuscript, 
including simulations, applications to spatial transcriptomics data, figure generation, and supplementary results.

## Previous working name and compatibility

**SpaRGS is the final method name.** The method was developed under the working name **sGSVA**. 
Some existing R functions, source files, saved objects, output filenames, and reproducibility metadata may
retain `sgsva` in their names. These historical names are preserved for backward compatibility and to maintain 
traceability to the original analyses; they refer to the SpaRGS method. 

## Project overview

SpaRGS maps pathway-associated expression contrasts across spatial locations without requiring a matched 
single-cell reference for scoring. The procedure estimates gene variability weights and directional subsets 
from **unsmoothed** expression, then aggregates ranks calculated from **spatially smoothed** expression. 
The workflow includes:

- Neighborhood-based Gaussian smoothing for rank calculation
- Correlation-based assignment of positive and negative pathway gene subsets
- Variability-weighted directional rank aggregation and within-pathway standardization
- Descriptive spatial organization using Moran's I and Geary's C
- Optional exploratory fixed-score randomization using global or within-hex-block permutations

Spatial autocorrelation and fixed-score randomization describe the organization of resulting maps; 
they do not establish calibrated inference for the complete expression-to-score pipeline. Biological annotations 
used in downstream evaluation are not required by the core scoring method.

The manuscript applications comprise:

- Breast cancer Visium
- Human dorsolateral prefrontal cortex (DLPFC) Visium
- Xenium FFPE breast cancer
- Simulation benchmarks and robustness analyses

## Repository structure

- `scripts/` — dataset- and analysis-specific workflows
- `R/` — reusable scoring, spatial-statistic, simulation, and evaluation functions
- `data/` — expected local input paths and the small frozen example input
- `results/` — selected machine-readable summaries, tables, figures, and analysis outputs
- `docs/` — manifests and supporting documentation

Public raw datasets and large intermediate analysis objects are not bundled with the lightweight repository unless explicitly indicated.

## Requirements

The manuscript workflows were developed using R 4.5.1. Dependencies vary by workflow; major packages are showing below. 
Additional packages and operating-system libraries may be needed for specific analyses.

### Install core analysis dependencies

```r
install.packages(c(
  "Matrix", "matrixStats", "spdep", "sf", "FNN", "msigdbr",
  "Seurat", "hexbin", "ggplot2", "patchwork", "dplyr",
  "tidyr", "readr", "MASS", "purrr", "tibble", "stringr",
  "pROC", "openxlsx", "jsonlite", "ggspatial", "ggpubr",
  "scales", "viridis", "reshape2", "gridExtra", "here", "fs"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
"GSVA", "org.Hs.eg.db", "spatialLIBD", "AUCell", "UCell",
"CelliD", "spacexr", "SpatialExperiment", "SummarizedExperiment",
"scater", "scran"
))

install.packages("remotes")

#install gsdensity from github
remotes::install_github("https://github.com/qingnanl/gsdensity")
```
Installing current package versions does not necessarily recreate the recorded manuscript environment. 
Consult the saved package-version records when reproducing the published analyses.

The exploratory breast Visium EnrichMap comparison additionally requires a Python environment named `enrichmap_env`. 
The recorded environment used Python 3.11.15, EnrichMap 0.2.3, NumPy 2.4.6, pandas 2.3.3, and AnnData 0.12.19. 
These versions are recorded in `results/breast_visium/reproducibility/breast_visium_enrichmap_python_package_versions.csv`.

This Python environment is required for the full manuscript pipeline but is not needed for the fixed-input minimal example.

The minimal fixed-input example needs only the packages used by its scoring and permutation scripts; 
it does not require downloading any real tissue dataset.


## Session information and random seeds

Session information and package versions are recorded in the relevant workflow outputs. Repository files 
such as `results/package_versions_sgsva.csv` may retain the previous working name for reproducibility. 
Stochastic analyses use script-specific seeds and settings; consult the corresponding analysis scripts 
and output metadata rather than assuming one seed for every workflow.

## Data access

### Public datasets

| Dataset | Source | Expected local path |
|---|---|---|
| Breast Visium | Zenodo record 4739739 | `data/breast_visium/` |
| DLPFC Visium | `spatialLIBD` package | Retrieved programmatically |
| Xenium FFPE breast cancer | 10x Genomics FFPE human breast cancer, 5K Human Pan-Tissue and Pathways panel plus 100 custom genes | `data/xenium/` |

Dataset access information and input requirements are described in Supplementary Table S16.

### Breast Visium

Place the required count matrices, spatial images and scale factors, tissue annotations, and breast single-cell 
reference files (GSE176078) under `data/breast_visium/`, using the paths expected by the corresponding scripts.

### Xenium FFPE breast cancer

Place the public Xenium expression matrix, cell centroid coordinates, and the gene-group marker table 
under `data/xenium/`, using the file structure expected by the analysis script.

### DLPFC

The DLPFC dataset is retrieved through `spatialLIBD` by the analysis scripts; no manual download into `data/` is required.

Raw datasets can be large and should be obtained from their original public sources before running full analyses.

## Key outputs

- Main and supplementary figures: `results/figures/`
- Machine-readable table sources: `results/tables/`
- Analysis-specific outputs: `results/simulation/`, `results/robustness/`, `results/breast_visium/`, `results/DLPFC/`, and `results/xenium/`
- Script-to-output manifest: `docs/script_to_output_manifest.xlsx`
- Minimal example output: `results/example/`

Supplementary Data 1–5 accompany the manuscript separately. The lightweight repository contains selected machine-readable summaries; 
it should not be assumed to include every full-resolution score matrix or permutation object.

## Manuscript parameter settings
| Analysis | Scoring graph | Bandwidth (σ) | Directionality threshold |
|---|---|---|---|
| Breast Visium | Delaunay | 150 µm | 0.1 |
| DLPFC | Delaunay | 151 µm | 0.2 |
| Xenium | KNN, k=8 | 24 µm | 0.1 |

These analyses used CV-weight clipping at 0.99, a minimum of 15 selected genes, a maximum of 25 genes per directional subset, 
and the `center_0.5` one-sided convention. Full spatial randomization used 1,000–10,000 permutations with seed 124. 
The minimal example uses 100 permutations. Analysis-specific settings and overrides are defined in the corresponding scripts.

## How to run

Run commands from the repository root after installing the required packages.

### Minimal reproducibility example

A lightweight fixed-input example is provided to test the core scoring and spatial-randomization workflow 
without downloading a tissue dataset. The expected frozen file is:

```text
data/example/minimal_example_input.rds
```

It contains simulated counts, spatial coordinates, representative KEGG Legacy pathways, latent truth fields, 
and simulation-design information. The input should **not be regenerated** during ordinary clean-clone testing; 
the one-time generator is for preparing the frozen input, not for running the example.

Install the minimal-example dependencies:

```bash
install.packages(c(
  "Matrix", "matrixStats", "FNN", "spdep", "dplyr"
))
```

From the repository root, run:

```bash
Rscript scripts/run_minimal_example.R
```

The script writes scoring, simulated-truth comparison, permutation smoke-test, and session-information 
outputs to `results/example/`. This is a functionality check with fixed inputs, **not** a reproduction 
of the manuscript's full simulation benchmark or spatial-inference settings. Some output filenames may still 
contain `sgsva` to preserve legacy interfaces and archived result references.

### Full manuscript pipeline

The full analysis workflows require public raw data and, for some steps, intermediate objects generated by preceding analyses. 
The repository documents the following pipeline entry point:

```bash
bash run_all.sh
```

Verify that this entry point and its prerequisites match the exact release being used. Individual workflows are also located under `scripts/`.


### Large figure-source files

Some manuscript figures require large intermediate objects, including Seurat objects, full pathway-score matrices, 
or permutation outputs. These are not included in the lightweight repository. When needed, regenerate them from 
the public input datasets using the corresponding full analysis scripts.

