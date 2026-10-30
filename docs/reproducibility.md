# SpaRGS reproducibility notes

This repository contains code and analysis workflows supporting the **SpaRGS** manuscript, including simulations, 
real-data applications, figures, tables, and selected machine-readable results. **sGSVA** was the previous working name; 
some existing functions, filenames, and archived outputs retain `sgsva` for backward compatibility and traceability.

## Entry points

The documented workflow entry points are:

- `bash run_all.sh`: full workflow, subject to the required public inputs and intermediate files
- `Rscript run_all_figures.R`: regenerate manuscript figures from the required analysis outputs
- `Rscript run_all_tables.R`: regenerate manuscript tables from the required analysis outputs
- `Rscript scripts/run_minimal_example.R`: run the fixed-input minimal example

The entry points and their prerequisites must be checked against the exact manuscript release; the minimal example 
does not reproduce the complete benchmark.

## Script-to-output manifest

`docs/script_to_output_manifest.xlsx` identifies the scripts and inputs used for the main and 
supplementary figures, tables, and Supplementary Data files.

## Random seeds and simulation generation

Random-number seeds are **analysis-specific**. There is no single seed that reproduces all SpaRGS simulations.

| Analysis step | Seed specification | Scope |
|---|---|---|
| Primary simulation benchmark | `base_seed = 10000L`; `scenario_seed = base_seed + (r - 1L) * nrow(grid) + g` | One seed for each replicate `r` and scenario row `g` |
| Spatial graph and hex-block robustness simulation | `base_seed = 2000L`; `set.seed(base_seed + i)` | One seed for each spatial-scale index `i` in `ell` |
| Hex-block robustness analysis | `set.seed(123 + p)` | One seed for each pathway index `p` |
| Frozen minimal-example input generation | `set.seed(1)` | Used only when creating the committed example input |
| Minimal-example permutation smoke test | `seed_permutation = 124` | Test-specific; does not use the manuscript permutation count |
| Full spatial randomization | Default `seed_permutation = 124` | Per invocation of the spatial-permutation workflow |
| Bootstrap resampling | `set.seed(1)` | The seeds are setting separately for simulation and Breast Visium analysis |
| Breast Visium, DLPFC, and Xenium analyses | `set.seed(1)` | Their seeds and stochastic steps are not determined by the simulation schedules above |

For the **primary simulation benchmark**, the seed is set immediately before each call to `simulate_counts()`. 
With 20 replicates and 18 scenario rows, the assigned seeds run from **10001 through 10360**. 
For the **spatial graph and hex-block robustness simulations**, the seed depends on the position of each spatial scale in `ell`:
For example, if `ell = c(50, 200, 400)`, the seeds are 2001, 2002, and 2003, respectively. The order and contents 
of `ell`, the construction and ordering of `grid`, and the number of replicates therefore form part of the 
reproducibility specification. Altering their order can change which seed is assigned to a condition.
The **Hex-block robustness analysis** assigned seeds from **124 through 309**.

## Dataset paths

| Dataset | Source | Expected local path |
|---|---|---|
| Breast Visium | Zenodo record 4739739 | `data/breast_visium/` |
| DLPFC Visium | `spatialLIBD` package | Retrieved programmatically |
| Xenium FFPE breast cancer | 10x Genomics FFPE human breast cancer dataset, 5K Human Pan-Tissue and Pathways panel plus 100 custom genes | `data/xenium/` |

Dataset access and expected input details are described in Supplementary Table S16. Large raw datasets are not redistributed in the lightweight repository.

## Minimal clean-installation test

The fixed-input example uses the committed file:

`data/example/minimal_example_input.rds`

The file was generated once from a reduced simulation design and is **not regenerated** in ordinary 
clean-clone testing. It contains the same counts, coordinates, gene sets, and simulated truth fields 
for every installation. 

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

This is a small installation and end-to-end workflow check, not a reproduction of the full simulation 
benchmark or manuscript-scale randomization. 

## Reproducibility limitations

Large public raw datasets must be obtained from their source repositories. Some full spatial-randomization 
analyses take several days on the reported hardware, and large intermediate objects are not included in the 
lightweight repository. Minor numerical differences can occur across software environments even when the input data and seeds agree.

## Session information

Session information and package versions are recorded in the workflow outputs, including `results/session_info.txt` 
and `results/package_versions_sgsva.csv` where present. The latter retains the previous working name to preserve 
traceability to archived results. Use the exact final release and recorded analysis environment when reproducing reported figures or tables.


