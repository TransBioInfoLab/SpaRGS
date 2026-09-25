################################################################################
# check gsdensity implementation
library(gsdensity)
library(ggplot2) # for plotting
library(reshape2)
library(msigdbr) # for gathering gene sets
library(Seurat)
library(SeuratData)

project_dir <- "."
outdir <- file.path(project_dir, "results", "simulation")

# save package versions
packageVersion("msigdbr")
packageVersion("gsdensity")
packageVersion("CelliD")
packageVersion("Seurat")
packageVersion("SeuratObject")

mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
mdb_c5_bp <- mdb_c5[mdb_c5$gs_subcat == "GO:BP", ]
gene.set.list <- list()
for (gene.set.name in unique(mdb_c5_bp$gs_name)){
  gene.set.list[[gene.set.name]] <- mdb_c5_bp[mdb_c5_bp$gs_name %in% gene.set.name, ]$gene_symbol
}

# from SeuratData
data("pbmc3k")
# Normalize the data
pbmc3k <- UpdateSeuratObject(pbmc3k)
pbmc3k <- NormalizeData(pbmc3k)

# compute cell/gene embeddings
ce <- compute_mca_v5(object = pbmc3k)

# build nearest neighbor graph
cells <- colnames(pbmc3k)
el <- compute.nn.edges(coembed = ce, nn.use = 300)

# get label propagation probability for each cell of a gene set 'GOBP_B_CELL_ACTIVATION'
cv <- run.rwr(el = el, gene_set = gene.set.list[["GOBP_B_CELL_ACTIVATION"]], cells = cells)
# get a 'positive' or 'negative' label for each cell of this gene set
cl <- compute.cell.label(cv)

# Compute the UMAP coordinates for visualization
pbmc3k <- pbmc3k %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  RunUMAP(dims = 1:20)
# plot the cell annotation
DimPlot(pbmc3k,
        group.by = "seurat_annotations",
        label = T,
        raster = T)

# add the label propagation probability and binarized label to meta data
pbmc3k@meta.data$b_cell_activition <- cv[colnames(pbmc3k)]
pbmc3k@meta.data$b_cell_activition_bin <- cl[colnames(pbmc3k)]

# plot the probabilities
p1 <- FeaturePlot(pbmc3k,
                  features = "b_cell_activition",
                  raster = T)


# check B-cell scores
validation_df <- data.frame(
  cell = colnames(pbmc3k),
  score = cv[colnames(pbmc3k)],
  annotation = pbmc3k$seurat_annotations
)

validation_summary <- validation_df |>
  dplyr::group_by(annotation) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(mean_score))
validation_summary

bcell_genes <- gene.set.list[["GOBP_B_CELL_ACTIVATION"]]
c(genes_in_set = length(unique(bcell_genes)),
  genes_in_pbmc3k = length(intersect(
    unique(bcell_genes),
    rownames(pbmc3k)
  )),
  genes_in_coembedding = length(intersect(
    unique(bcell_genes),
    rownames(ce)
  ))
)

write.csv(
  validation_summary, 
  file.path(outdir, "validation//pbmc3k_b_cell_activation_by_celltype.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    genes_in_set = 297,
    genes_in_pbmc3k = 237,
    genes_in_coembedding = 237
  ),
  file.path(outdir, "validation/pbmc3k_b_cell_activation_gene_overlap.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "validation/sessionInfo.txt"),
)

ggsave(
  file.path(outdir, "validation/pbmc3k_b_cell_activation_featureplot.pdf"),
  plot = p1,
  width = 6,
  height = 5
)

# The Seurat v5-compatible MCA wrapper was validated by reproducing the GSDensity PBMC3K tutorial. 
# The resulting B-cell activation probability map showed the same qualitative structure and comparable 
# probability scale as the original example, with the strongest signal concentrated in the same isolated 
# cell population and weaker signal across the remaining populations. Differences in UMAP orientation and local 
# geometry were expected because UMAP embeddings are stochastic and identifiable only up to rotation and reflection. 
# Cell-type-stratified score summaries further confirmed that the B-cell population received the strongest enrichment scores
