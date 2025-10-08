suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(reticulate)
  library(sceasy)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(anndata)
})

csv_path <- "config/proyects_metadata_samples_fastqs.csv"
root_counts <- "/datos/sensence/emilio/liver_sc/fastq_prepros/output/counts"
out_root <- "results/macrophages"
layer <- "filtered"

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# ---- READ METADATA ----
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- unique(df[, c("proyect", "ident_sample")])

proj_col <- df$proyect
samp_col <- df$ident_sample

# ---- HELPERS ----
mtx_dir_for <- function(project, sample, layer) {
  file.path(
    root_counts, project,
    paste0(sample, "_matrix_Solo.out"), "Gene", layer
  )
}

# ---- MERGING ----
obj_list <- list()
for (i in seq_len(nrow(df))) {
  project <- proj_col[i]
  sample <- samp_col[i]
  mtx_dir <- mtx_dir_for(project, sample, layer)
  if (!dir.exists(mtx_dir)) {
    message(sprintf("SKIP: %s/%s -- missing %s", project, sample, mtx_dir))
    next
  }
  message(sprintf("[Loading] %s / %s layer=%s", project, sample, layer))
  # see if it is gzipped as read10x needs it
  is_gz <- function(path) {
    if (!file.exists(path)) stop("No such file: ", path)
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    sig <- readBin(con, "raw", 2L)
    length(sig) == 2L && as.integer(sig[1]) == 0x1f && as.integer(sig[2]) == 0x8b
  }
  # Compress in-place with system gzip if needed; returns gz path
  gzip_if_needed <- function(path) {
    stopifnot(file.exists(path), !dir.exists(path))
    if (is_gz(path)) {
      if (!grepl("\\.gz$", path)) { # gz content but no .gz suffix → rename
        dest <- paste0(path, ".gz")
        file.rename(path, dest)
        return(dest)
      }
      return(path) # already gz and named *.gz
    }
    # plain file; if name ends with .gz but content isn't gz, fix name first
    if (grepl("\\.gz$", path)) {
      base <- sub("\\.gz$", "", path)
      file.rename(path, base)
      path <- base
    }
    system2("gzip", c("-f", shQuote(path)), stdout = NULL, stderr = NULL)
    paste0(path, ".gz")
  }
  gz_paths <- sapply(
    c("matrix.mtx", "barcodes.tsv", "features.tsv"),
    function(f) {
      f <- file.path(mtx_dir, f)
      if (file.exists(f)) gzip_if_needed(f)
    },
    USE.NAMES = TRUE
  )
  # Load it to seurat since regex wont work now with Read10x function
  pick <- function(d, a, b) if (file.exists(file.path(d, a))) file.path(d, a) else file.path(d, b)
  mtx <- pick(mtx_dir, "matrix.mtx.gz", "matrix.mtx")
  bc <- pick(mtx_dir, "barcodes.tsv.gz", "barcodes.tsv")
  feat <- {
    f <- pick(mtx_dir, "features.tsv.gz", "features.tsv")
    if (!file.exists(f)) f <- pick(mtx_dir, "genes.tsv.gz", "genes.tsv")
    f
  }
  # Count columns in the first line of features to choose the right name column
  feat_cols <- length(strsplit(readLines(feat, n = 1), "\t", fixed = TRUE)[[1]])
  gene_col <- if (feat_cols >= 2) 2 else 1
  # load seurat obj
  counts <- ReadMtx(
    mtx = mtx,
    cells = bc,
    features = feat,
    feature.column = gene_col
  )
  obj <- CreateSeuratObject(
    counts = counts,
    project = paste(project, sample, sep = "_")
  )
  obj$project <- project
  obj$sample <- sample
  obj_list[[paste(project, sample, sep = "_")]] <- obj
}

hca_liver_raw <- Reduce(f = merge, obj_list)
# correct the merge seurat proyect duplucation naming on the layers
hca_liver_raw <- JoinLayers(hca_liver_raw, assay = "RNA")
hca_liver_raw[["RNA"]] <- split(
  hca_liver_raw[["RNA"]],
  hca_liver_raw$orig.ident
)

# ---- QC ----
hca_liver_raw <- NormalizeData(hca_liver_raw, verbose = FALSE)
hca_liver_raw <- JoinLayers(hca_liver_raw, assay = "RNA")
hca_liver_raw <- PercentageFeatureSet(hca_liver_raw,
  pattern = "^MT-",
  col.name = "percent.mt"
)
hca_liver_raw <- subset(hca_liver_raw, subset = percent.mt < 30) # nolint
hca_liver_raw[["RNA"]] <- split(
  hca_liver_raw[["RNA"]],
  hca_liver_raw$orig.ident
)
hca_liver_raw <- FindVariableFeatures(hca_liver_raw)
hca_liver_raw <- ScaleData(hca_liver_raw)
hca_liver_raw <- RunPCA(hca_liver_raw, npcs = 50)

# ---- SCVI integration ----
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
hca_liver_scvi <- JoinLayers(hca_liver_raw, assay = "RNA")
sce <- as.SingleCellExperiment(hca_liver_scvi, assay = "RNA")
out_sce <- file.path(out_root, "hca_liver_scvi.h5ad")
writeH5AD(sce, out_sce)
adata <- anndata::read_h5ad(file.path(out_root, "hca_liver_scvi.h5ad"))
print(adata) # Note generally in Python, dataset conventions are obs x var
# run setup_anndata
scvi$model$SCVI$setup_anndata(adata, batch_key = "sample")
# create the model
model <- scvi$model$SCVI(adata)
# train the model
model$train()
latent <- model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) <- colnames(hca_liver_scvi)
hca_liver_scvi[["scvi"]] <- CreateDimReducObject(
  embeddings = latent,
  key = "scvi_",
  assay = DefaultAssay(hca_liver_scvi)
)
hca_liver_scvi <- FindNeighbors(hca_liver_scvi, dims = 1:10, reduction = "scvi")
hca_liver_scvi <- FindClusters(hca_liver_scvi, resolution = 1)

hca_liver_scvi <- RunUMAP(hca_liver_scvi,
  dims = 1:10,
  reduction = "scvi",
  n.components = 2,
  reduction.name = "scvi_umap"
)
DimPlot(hca_liver_scvi, reduction = "umap", pt.size = 3)
p1 <- DimPlot(hca_liver_scvi, reduction = "scvi_umap", group.by = "orig.ident", pt.size = 2)
p1
saveRDS(hca_liver_scvi,
  file = file.path(out_root, "hca_liver_scvi.rds")
)

# ---- Harmony integration ----
hca_liver_har[["RNA"]] <- split(
  hca_liver_raw[["RNA"]],
  hca_liver_raw$orig.ident
)
hca_liver_har <- RunPCA(hca_liver_raw, npcs = 50)
hca_liver_har <- IntegrateLayers(hca_liver_har,
  HarmonyIntegration,
  new.reduction = "harmony"
)
hca_liver_har <- FindNeighbors(hca_liver_har,
  reduction = "harmony",
  dims = 1:30
)
hca_liver_har <- FindClusters(hca_liver_har,
  resolution = 0.5,
  cluster.name = "harmony_clusters"
)
hca_liver_har <- RunUMAP(hca_liver_har,
  reduction = "harmony",
  dims = 1:50, reduction.name = "harmony_umap"
)
DimPlot(hca_liver_har, reduction = "harmony_umap", pt.size = 3)
saveRDS(hca_liver_har,
  file = file.path(out_root, "hca_liver_har.rds")
)

# ---- Macrphage exploring ----
macrophages <- subset(hca_liver_raw, cells = WhichCells(
  hca_liver_raw,
  expression =
    (CD68 > 0.5 |
      ADGRE1 > 0.5 |
      ITGAM > 0.5 |
      CSF1R > 0.5 |
      MERTK > 0.5 |
      FCGR1A > 0.5 |
      MARCO > 0.5) &
      PTPRC > 0.5,
  slot = "data"
))

macrophages <- NormalizeData(macrophages)
macrophages <- JoinLayers(macrophages, assay = "RNA") # nolint
macrophages <- PercentageFeatureSet(macrophages, pattern = "^MT-", col.name = "percent.mt")
macrophages <- subset(macrophages, subset = percent.mt < 30) # nolint
macrophages[["RNA"]] <- split(macrophages[["RNA"]], macrophages$orig.ident)
macrophages <- FindVariableFeatures(macrophages)
macrophages <- ScaleData(macrophages)
macrophages <- RunPCA(macrophages, npcs = 50)
sc <- import("scanpy", convert = FALSE)
