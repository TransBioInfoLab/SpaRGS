# Reproducibility notes

This repository contains code and workflows used to reproduce the analyses, figures, and tables for the sGSVA manuscript.

## Entry points

- `bash run_all.sh`: full workflow
- `Rscript run_all_figures.R`: regenerate all manuscript figures from processed outputs
- `Rscript run_all_tables.R`: regenerate all manuscript tables from processed outputs
- `Rscript scripts/run_minimal_example.R`: run a small example

## Script-to-output manifest

The file `docs/script_to_output_manifest.xlsx` maps each main figure, supplementary figure, main table, supplementary table, 
and supplementary Data files to the script and input files used to generate it.

## Random seeds

Fixed seeds are used for stochastic steps:

| Analysis step | Seed |
|---|---:|
| Reference minimal-example input generation | 1 |
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
| Xenium FFPE breast cancer | 10x Genomics public Xenium FFPE human breast cancer dataset with 5k Human Pan-Tissue and Pathways panel plus 100 Custom Genes | `data/xenium/` |

Exact dataset accessions and vendor download links are provided in Supplementary Table S16.

## Minimal clean-installation test

The minimal example uses the committed frozen input:

`data/example/minimal_example_input.rds`

This file was generated once from a reduced version of the simulation framework
and is not regenerated during ordinary clean-installation testing. Freezing the
input ensures that independent systems begin from the same counts, spatial
coordinates, pathway definitions, and simulated truth fields.

From the repository root, run:

```bash
Rscript scripts/run_minimal_example.R
```

## Reproducibility limitations

Large raw datasets are not redistributed. Users must obtain the public datasets
from their original sources and place them in the expected local directories.
Some full permutation analyses are computationally intensive and may require
multiple days depending on hardware and permutation settings.

The minimal example is intended as a lightweight installation and workflow
check rather than a reproduction of the full manuscript benchmark. It uses a
committed frozen input so that independent installations begin from identical
data. Small numerical differences across operating systems or numerical
libraries may still occur and should be interpreted within numerical tolerance.

## Session information

Session information is saved in `results/session_info.txt`. A curated list of packages used in the sGSVA workflow is saved in `results/package_versions_sgsva.csv`.
