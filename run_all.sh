#!/bin/bash
set -eu

echo "Running a minimal example..."
Rscript scripts/run_example_visium.R
Rscript scripts/run_example_xenium.R

echo "Running simulation analyses..."
Rscript scripts/Simulation/simulation_analysis.R
Rscript scripts/Simulation/simulation_GSDesnity.R

echo "Running robustness analyses..."
Rscript scripts/Robustness/robustness_analysis.R

echo "Running Breast Visium analyses..."
Rscript scripts/BreastVisium/breast_visium_with_gsdensity.R
Rscript scripts/BreastVisium/breast_celltype_anchoring_all_methods.R
Rscript scripts/BreastVisium/breast_visium_enrichmap.R

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