#!/bin/bash
set -eu

echo "Running a minimal example..."
Rscript scripts/run_minimal_example.R

echo "Running simulation analyses..."
Rscript scripts/Simulation/simulation_single_sample_geneset.R
Rscript scripts/Simulation/simulation_with_gsdensity.R

echo "Running robustness analyses..."
Rscript scripts/Robustness/robustness_analysis.R

echo "Running Breast Visium analyses..."
Rscript scripts/BreastVisium/breast_visium_with_gsdensity.R
Rscript scripts/BreastVisium/breast_visium_celltype_anchoring_all_methods.R
Rscript scripts/BreastVisium/01_breast_visium_export_enrichmap_input.R

conda activate enrichmap_env
python scripts/BreastVisium/02_run_enrichmap_breast_visium.py

Rscript scripts/BreastVisium/03_breast_visium_evaluate_enrichmap.R

echo "Running DLPFC analyses..."
Rscript scripts/DLPFC/DLPFC_analysis.R

echo "Running Xenium analyses..."
Rscript scripts/Xenium/Xenium_analysis.R

echo "Generating all figures..."
Rscript run_all_figures.R

echo "Generating all tables..."
Rscript run_all_tables.R

echo "Saving session information..."
Rscript scripts/session_info.R

echo "All analyses completed successfully."