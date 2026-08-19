# Reproducibility notes

This repository contains code and workflows used to reproduce the analyses, figures, and tables for the sGSVA manuscript.

## Entry points

- `bash run_all.sh`: full workflow
- `Rscript run_all_figures.R`: regenerate all manuscript figures from processed outputs
- `Rscript run_all_tables.R`: regenerate all manuscript tables from processed outputs
- `Rscript scripts/run_example_visium.R`: run a small Breast Visium example
- `Rscript scripts/run_example_xenium.R`: run a small Xenium FFPE example

## Script-to-output manifest

The file `docs/script_to_output_manifest.xlsx` maps each main figure, supplementary figure, main table, supplementary table, and supplementary Data files to the script and input files used to generate it.

## Random seeds

Fixed seeds are used for stochastic steps:

| Analysis step | Seed |
|---|---:|
| Minimal example sampling | 1 |
| Simulation benchmark | 1 |
| Robustness simulations | 1 |
| Breast Visium analysis | 1 |
| DLPFC analysis | 1 |
| Xenium analysis | 1 |
| Spatial permutations | 124 |

## Dataset paths

| Dataset | Source | Expected local path |
|---|---|---|
| Breast Visium | Zenodo accession 4739739 | `data/breast_visium/` |
| DLPFC Visium | `spatialLIBD` package | downloaded programmatically |
| Xenium FFPE breast cancer | 10x Genomics public Xenium 5k Pan-Tissue and Pathways dataset | `data/xenium/` |
Exact dataset accessions and vendor download links are provided in Supplementary Table S16.

## Reproducibility limitations

Large raw datasets are not redistributed. Users must download the public datasets and place them in the expected folders. Some long-running permutation analyses may require multiple days depending on hardware and permutation settings.

## Session information

Session information is saved in `results/session_info.txt`. A curated list of packages used in the sGSVA workflow is saved in `results/package_versions_sgsva.csv`.