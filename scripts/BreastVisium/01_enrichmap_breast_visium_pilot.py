#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import json
import time
import platform
import importlib.metadata as md
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import enrichmap as em


# In[2]:


# ------------------------------------------------------------
# User configuration
# ------------------------------------------------------------

from pathlib import Path
import os

# Option 1: assume the notebook is run from the project root
project_dir = Path.cwd()

# If the notebook is inside scripts/BreastVisium, move up to project root
if project_dir.name == "BreastVisium":
    project_dir = project_dir.parents[1]
elif project_dir.name == "scripts":
    project_dir = project_dir.parent

pilot_dir = project_dir / "results" / "breast_visium" / "enrichmap_pilot"
prefix = "breast_visium"
n_neighbors = 6

expr_file = pilot_dir / f"{prefix}_expr_spots_by_genes.csv"
coords_file = pilot_dir / f"{prefix}_coords.csv"
geneset_file = pilot_dir / f"{prefix}_gene_sets.json"

print("Project directory:", project_dir)
print("Pilot directory:", pilot_dir)
print("Expression file exists:", expr_file.exists())
print("Coordinates file exists:", coords_file.exists())
print("Gene-set file exists:", geneset_file.exists())


# In[3]:


# ------------------------------------------------------------
# Load exported breast Visium input
# ------------------------------------------------------------

expr = pd.read_csv(expr_file).set_index("spot_id")
coords = pd.read_csv(coords_file).set_index("spot_id")

# Make sure coordinates are in the same spot order as expression
coords = coords.loc[expr.index]

with open(geneset_file, "r") as f:
    gene_sets_raw = json.load(f)

# Keep only genes present in the expression matrix
gene_sets = {
    k: [g for g in v if g in expr.columns]
    for k, v in gene_sets_raw.items()
}

gene_set_coverage = pd.DataFrame({
    "pathway": list(gene_sets_raw.keys()),
    "n_genes_original": [len(v) for v in gene_sets_raw.values()],
    "n_genes_in_expr": [len(gene_sets[k]) for k in gene_sets_raw.keys()]
})

gene_set_coverage["passed_min_genes"] = gene_set_coverage["n_genes_in_expr"] >= 2

gene_set_coverage.to_csv(
    pilot_dir / f"{prefix}_enrichmap_gene_set_coverage_from_python.csv",
    index=False
)

print("Expression:", expr.shape)
print("Coordinates:", coords.shape)
print("Gene sets:", len(gene_sets))
display(gene_set_coverage.head())


# In[5]:


# ------------------------------------------------------------
# Create AnnData object
# ------------------------------------------------------------

adata = ad.AnnData(
    X=expr.values,
    obs=pd.DataFrame(index=expr.index),
    var=pd.DataFrame(index=expr.columns)
)

adata.obsm["spatial"] = coords[["x", "y"]].values

adata


# In[6]:


# ------------------------------------------------------------
# Run EnrichMap in two modes
# ------------------------------------------------------------

summary_rows = []

enrichmap_modes = {
    "default_corrected": True,
    "smooth_only": False
}

for mode_name, correct_spatial in enrichmap_modes.items():

    print(f"Running EnrichMap mode: {mode_name}")

    adata_mode = adata.copy()

    # Prefix pathway names so EnrichMap output columns are unique by mode
    gs_mode = {
        f"{mode_name}__{k}": v
        for k, v in gene_sets.items()
    }

    t0 = time.time()

    em.tl.score(
        adata=adata_mode,
        gene_set=gs_mode,
        spatial_key="spatial",
        n_neighbors=n_neighbors,
        smoothing=True,
        correct_spatial_covariates=correct_spatial,
        batch_key=None
    )

    runtime = time.time() - t0

    score_cols = [
        c for c in adata_mode.obs.columns
        if c.endswith("_score")
    ]

    print(mode_name, "score columns:", len(score_cols))
    print(mode_name, "runtime seconds:", runtime)

    score_df = adata_mode.obs[score_cols].copy()

    rename_map = {}
    for c in score_cols:
        name = c.replace("_score", "")
        name = name.replace(f"{mode_name}__", "")
        rename_map[c] = name

    score_df = score_df.rename(columns=rename_map)

    # Convert to pathway x spot
    score_mat = score_df.T

    score_outfile = pilot_dir / f"{prefix}_enrichmap_{mode_name}_all_scores.csv"
    score_mat.to_csv(score_outfile)

    summary_rows.append({
        "dataset": prefix,
        "mode": mode_name,
        "runtime_seconds": runtime,
        "n_spots": adata_mode.n_obs,
        "n_genes": adata_mode.n_vars,
        "n_gene_sets_input": len(gene_sets),
        "n_score_columns": len(score_cols),
        "n_pathways_scored": score_mat.shape[0],
        "n_neighbors": n_neighbors,
        "smoothing": True,
        "correct_spatial_covariates": correct_spatial,
        "score_output": str(score_outfile)
    })

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    pilot_dir / f"{prefix}_enrichmap_all_runtime_summary.csv",
    index=False
)

display(summary)


# In[7]:


# ------------------------------------------------------------
# Save Python environment information
# ------------------------------------------------------------

def safe_version(pkg):
    try:
        return md.version(pkg)
    except md.PackageNotFoundError:
        return None

versions = {
    "python_executable": sys.executable,
    "python": platform.python_version(),
    "numpy": safe_version("numpy"),
    "pandas": safe_version("pandas"),
    "anndata": safe_version("anndata"),
    "enrichmap": safe_version("enrichmap")
}

version_df = pd.DataFrame(
    [{"package": k, "version": v} for k, v in versions.items()]
)

version_df.to_csv(
    pilot_dir / f"{prefix}_enrichmap_python_package_versions.csv",
    index=False
)

display(version_df)

