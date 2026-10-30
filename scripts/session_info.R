# scripts/session_info.R
# Save R session information and package versions for the SpaRGS analysis.

session_file <- file.path("results", "session_info.txt")
package_file <- file.path("results", "package_versions_sgsva.csv")

# Packages used by the SpaRGS manuscript workflow.
sgsva_packages <- c(
  "Seurat",
  "GSVA",
  "AUCell",
  "UCell",
  "msigdbr",
  "spdep",
  "sf",
  "FNN",
  "hexbin",
  "Matrix",
  "matrixStats",
  "ggplot2",
  "patchwork",
  "pheatmap",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "purrr",
  "stringr",
  "data.table",
  "spatialLIBD",
  "SpatialExperiment",
  "org.Hs.eg.db",
  "scater",
  "spacexr",
  "pROC",
  "openxlsx",
  "jsonlite"
)

# Save sessionInfo for the current environment.
writeLines(
  c(
    paste0("Generated: ", Sys.time()),
    paste0("Working directory: ", normalizePath(getwd())),
    "",
    capture.output(sessionInfo())
  ),
  con = session_file
)

# Save versions only for packages relevant to the SpaRGS workflow.
installed <- installed.packages()

pkg_versions <- data.frame(
  Package = sgsva_packages,
  Version = vapply(
    sgsva_packages,
    function(pkg) {
      if (pkg %in% rownames(installed)) {
        as.character(packageVersion(pkg))
      } else {
        NA_character_
      }
    },
    character(1)
  ),
  Installed = sgsva_packages %in% rownames(installed),
  stringsAsFactors = FALSE
)

write.csv(pkg_versions, package_file, row.names = FALSE)

message("Wrote session information to: ", session_file)
message("Wrote SpaRGS package versions to: ", package_file)

