# Third party notices

This file records identified third-party source adaptations and the external pathway-scoring dependencies used by the SpaRGS analysis workflows. Source and installed-package metadata were checked on 2026-11-02.

## License scope

The root [LICENSE](LICENSE) contains the MIT license for original, independently authored SpaRGS code. It does not replace the licenses of third-party code or grant MIT-only rights to the GSDensity-derived functions identified below.

The GSDensity adaptations retain GNU General Public License version 3 terms (`GPL-3.0-only`, declared as `GPL-3` in the upstream R package). A complete, unmodified copy of that license is included at [LICENSES/GPL-3.0-only.txt](LICENSES/GPL-3.0-only.txt).

The GPL terms govern distribution of covered adaptations and any covered combined work. Identifying particular functions here does not limit the scope of those terms or establish that placing code in a separate file avoids them. Preserve applicable upstream notices, identify modifications, include the license, and provide the corresponding source required by the GPL when distributing a covered work. This provenance inventory is not a repository-wide license-clearance certificate.

## Adapted GSDensity source

### Upstream reference

- Project: GSDensity (`gsdensity`).
- Upstream author identified in the consulted repository's `DESCRIPTION`: Qingnan Liang.
- Repository consulted: <https://github.com/KChen-lab/gsdensity>.
- Source: `compute.mca()` in [`R/gsdensity_functions.R`](https://github.com/KChen-lab/gsdensity/blob/master/R/gsdensity_functions.R).
- License declaration: [`DESCRIPTION`](https://github.com/KChen-lab/gsdensity/blob/master/DESCRIPTION), `License: GPL-3`.
- Installed reference used for comparison: `gsdensity` 0.1.0, remote repository `KChen-lab/gsdensity`, recorded revision `d6c938e6d53e195f87d10cfdef88d6e766e0c6ab`.

The recorded installed revision identifies the reference available for this review; it does not establish the date or exact revision originally consulted when the adaptations were written. Upstream authorship and rights are retained. No upstream copyright year has been invented for a function that had no separate copyright-year notice in the inspected source.

### Affected functions

Both of these files define a function named `compute_mca_v5()` adapted from the upstream `compute.mca()`:

- [`scripts/Simulation/simulation_with_gsdensity.R`](scripts/Simulation/simulation_with_gsdensity.R).
- [`scripts/BreastVisium/breast_visium_with_gsdensity.R`](scripts/BreastVisium/breast_visium_with_gsdensity.R).

The adaptations preserve the upstream sequence of intersecting genes, normalizing expression, calculating MCA, and combining gene and cell coordinates. Changes relative to the inspected upstream function include:

1. Explicit assay selection and Seurat v5 layer access through `SeuratObject::LayerData()`.
2. Subsetting expression before conversion to a dense matrix.
3. Calling the matrix interface of `CelliD::RunMCA()` and extracting `featuresCoordinates` and `cellsCoordinates`.
4. Preserving matrix dimensions during coordinate selection with `drop = FALSE`.
5. In the simulation variant, exposing the normalization scale factor as an argument, with an existing default of `1e6`, and passing the selected assay explicitly to normalization. The breast Visium variant uses the default `NormalizeData()` normalization settings.

Attribution and modification notices were added by the SpaRGS contributors on 2026-11-02. The compatibility adaptations already existed before that documentation change; their original modification dates were not recorded. Adding these notices did not change either function's executable statements, parameters, or output handling.

The adapted functions are distributed under GNU GPL version 3, without any warranty, including the implied warranties of merchantability or fitness for a particular purpose. See the full license for the complete terms. The MIT license for original SpaRGS code does not relicense these adaptations.

## Academic attribution and other resources

Cite the scoring methods and software used in the manuscript and downstream work. For R packages, `citation("packageName")` returns the installed package's suggested citation. EnrichMap's repository provides its method citation. Scientific citations complement the source-code notices and do not replace license obligations.

Public tissue datasets, annotations, images, and MSigDB/KEGG gene-set memberships have separate attribution and reuse terms. The MIT software license and the GPL copy included here do not relicense those resources. Retain their source acknowledgements and applicable notices in the data archive.
