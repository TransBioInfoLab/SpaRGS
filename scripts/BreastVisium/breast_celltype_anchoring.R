# cell-type anchoring for breast cancer Visium

suppressPackageStartupMessages({
  library(annotate); library(frma); library(affy)
  library(Matrix); library(spacexr); library(Seurat);
  library(dplyr); library(ggplot2)
})

source("R/sgsva_score.R")

set.seed(1)

# ---------------------------
# paths
# ---------------------------
project_dir <- "."
data_dir    <- file.path(project_dir, "data", "breast_visium")
results_dir <- file.path(project_dir, "results", "breast_visium")

# sgsva file
sgsva_file <- file.path(results_dir, "sgsva_scores.csv")

# TNBC single-cell reference
sample <- "CID4495"
# Define file paths
metadata_file <- file.path(data_dir, "GSE176078/extracted/", sample, "metadata.csv")
matrix_file <- file.path(data_dir, "GSE176078/extracted/", sample, "count_matrix_sparse.mtx")
genes_file <- file.path(data_dir, "GSE176078/extracted/", sample, "count_matrix_genes.tsv")
barcodes_file <- file.path(data_dir, "GSE176078/extracted/", sample, "count_matrix_barcodes.tsv")


# ---------------------------
# load breast object + scores
# ---------------------------
bc_obj <- readRDS(file.path(results_dir, "bc_obj.rds"))

sgsva_scores <- read.csv(sgsva_file, row.names = 1, check.names = FALSE)
z_sgsva <- z_normalization(as.matrix(sgsva_scores))
t_sgsva <- as.data.frame(t(z_sgsva))
t_sgsva$barcode <- rownames(t_sgsva)


# ---------------------------
# load single-cell reference
# ---------------------------
# Read data
expression_matrix <- readMM(matrix_file)
barcodes <- read.table(barcodes_file, header = FALSE)$V1
genes <- read.table(genes_file, header = FALSE)
metadata <- read.csv(metadata_file, row.names = 1)

# Assign row/column names
colnames(expression_matrix) <- barcodes
rownames(expression_matrix) <- genes$V1

# Create Seurat object
sc_obj <- CreateSeuratObject(counts = expression_matrix, project = sample)
sc_obj <- AddMetaData(sc_obj, metadata)

# Process data
sc_obj <- SCTransform(sc_obj, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:20)
Idents(sc_obj) <- sc_obj@meta.data$celltype_major
DimPlot(sc_obj, reduction = "umap", label = TRUE)


# ---------------------------
# build RCTD reference
# ---------------------------
# Create a SpaceXr reference object
sc_counts <- GetAssayData(sc_obj, layer = "counts")
cell_types <- sc_obj@meta.data$celltype_major
names(cell_types) <- barcodes
cell_types <- as.factor(cell_types)
reference <- Reference(sc_counts, cell_types)

# create a SpaceXr spatial object
spatial_counts <- GetAssayData(bc_obj, assay = "Spatial", layer="counts")
coords <- GetTissueCoordinates(bc_obj) %>% as.data.frame()
coords$cell <- NULL
puck <- SpatialRNA(coords, spatial_counts)

# Create RCTD object
rctd <- create.RCTD(puck, reference, max_cores = 4)
rctd <- run.RCTD(rctd)
cell_type_proportions <- rctd@results$weights %>% as.data.frame()

# merge with sgsva scores
cell_type_proportions$barcode <- rownames(cell_type_proportions)
merge_dat <- inner_join(cell_type_proportions, t_sgsva, by="barcode")
write.csv(
  merge_dat,
  file.path(results_dir, "anchoring_sgsva_scores_merged.csv"),
  row.names = FALSE
)


# ---------------------------
# pathway-celltype correlation table
# ---------------------------
pathways <- colnames(t_sgsva)[!colnames(t_sgsva) %in% "barcode"]
celltypes <- colnames(cell_type_proportions)[!colnames(cell_type_proportions) %in% "barcode"]

cor_table <- tidyr::crossing(pathway = pathways, cell_type = celltypes) %>%
  rowwise() %>%
  mutate(cor = cor(merge_dat[[pathway]], merge_dat[[cell_type]], method = "spearman"),
         pval = cor.test(merge_dat[[pathway]], merge_dat[[cell_type]], method = "spearman")$p.value) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  arrange(desc(abs(cor)))
write.csv(cor_table, file.path(results_dir, "correlation_score_vs_celltype.csv"))

