reticulate::py_discover_config()

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
mdmeta_path <- "config/liver_data.csv"
root_counts <- "/datos/sensence/emilio/liver_sc/fastq_prepros/output/counts"
out_root <- "results/macrophages"
layer <- "filtered"

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# ---- READ METADATA ----
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
df <- unique(df[, c("proyect", "ident_sample")])

proj_col <- df$proyect
samp_col <- df$ident_sample

md_meta <- read.csv(mdmeta_path, header = TRUE)

# ---- HELPERS ----
mtx_dir_for <- function(project, sample, layer) {
  file.path(
    root_counts, project,
    paste0(sample, "_matrix_Solo.out"), "Gene", layer
  )
}
md_meta <- md_meta |>
  dplyr::select(-source_metadata, -cells_estimate, -source_data) |>
  dplyr::rename(project = proyect,
                sample = ident)
#grab the middle value of the intervals
md_meta <- md_meta %>%
  mutate(
    age = if_else(
      # condition: contains a dash, i.e., it's a range
      str_detect(age, "-"),
      
      # TRUE case → compute midpoint
      {
        # extract the two numbers from "min-max"
        parts <- str_split(age, "-", simplify = TRUE)
        mid <- (as.numeric(parts[,1]) + as.numeric(parts[,2])) / 2
        as.character(mid)
      },
      
      # FALSE case → leave as is
      age
    ),
    # finally convert to numeric
    age = as.numeric(age)
  )
md_meta <- md_meta %>%
  mutate(
    aging = if_else(age >= 50, "aged", "young")
  )
md_meta <- md_meta %>%
  filter(project != "chan_zuckerberg")
#classify it as young or aged on cutoff 50 years

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
hca_liver_scvi <- FindClusters(hca_liver_scvi, resolution = 0.5)

hca_liver_scvi <- RunUMAP(hca_liver_scvi,
  dims = 1:10,
  reduction = "scvi",
  n.components = 2,
  reduction.name = "scvi_umap"
)
DimPlot(hca_liver_scvi, reduction = "scvi_umap", pt.size = 3)
p1 <- DimPlot(hca_liver_scvi, reduction = "scvi_umap", group.by = "orig.ident", pt.size = 2)
p1
saveRDS(hca_liver_scvi,
  file = file.path(out_root, "hca_liver_scvi.rds")
)

# ---- Harmony integration ----
hca_liver_har <- hca_liver_raw
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

# ---- Macrphage exploring har ----
FeaturePlot(
  hca_liver_har, #or hca_liver_har
  features = c("CD68", "ADGRE1", "ITGAM", "CSF1R", "MERTK", "FCGR1A", "MARCO", "PTPRC"),
  reduction = "scvi_umap",
  min.cutoff = "q05",
  max.cutoff = "q95",
  order = TRUE, # plot high-expressers on top
  slot = "data" # use normalized log1p data for plotting
)
p <- DimPlot(
  hca_liver_har, reduction = "harmony_umap",
  label = TRUE, repel = TRUE, label.size = 5,
  raster = TRUE
) + NoLegend()

# ---- Macrophage obj cleaning ----
macrophages <- subset(
  hca_liver_har,
  subset = harmony_clusters %in% c(7, 8, 10, 17, 23)
)

#reprocess it
DefaultAssay(macrophages) <- "RNA"
macrophages <- NormalizeData(macrophages)
macrophages <- FindVariableFeatures(macrophages, nfeatures = 3000)
macrophages <- ScaleData(macrophages, features = VariableFeatures(macrophages))
macrophages <- RunPCA(macrophages, features = VariableFeatures(macrophages))
macrophages <- IntegrateLayers(macrophages,
                               HarmonyIntegration,
                               new.reduction = "harmony")
macrophages <- FindNeighbors(macrophages,
                             reduction = "harmony",
                             dims = 1:30)
macrophages <- FindClusters(macrophages, resolution = 0.3)
macrophages <- RunUMAP(macrophages,
                       reduction = "harmony",
                       reduction.name = "harmony_umap",
                       dims = 1:30)
try(httpgd::hgd_close(all = TRUE), silent = TRUE)
if (!is.null(dev.list())) graphics.off()
httpgd::hgd(port = 0)
httpgd::hgd_browse()
DimPlot(
  macrophages, reduction = "harmony_umap",
  label = TRUE, repel = TRUE, label.size = 5,
  raster = TRUE
) + NoLegend()
FeaturePlot(macrophages,
            features = c("TIMD4","MARCO","VSIG4","CD5L","TREM2","CD9","SPP1"),
            reduction = "harmony_umap", min.cutoff = "q05", max.cutoff = "q95", order = TRUE)

Idents(macrophages) <- "seurat_clusters"
macrophages <- JoinLayers(macrophages, assay = "RNA")
library(presto)
markers <- FindAllMarkers(macrophages, only.pos = TRUE, logfc.threshold = 0.25)
head(markers[order(markers$avg_log2FC, decreasing = TRUE), ], 20)

FeaturePlot(macrophages,
            features = c("CDKN1A", "CDKN2A"),
            reduction = "harmony_umap",
            min.cutoff = "q05",
            max.cutoff = "q95",
            order = TRUE)

#missing a solid anotation for marophage subtype

macrophages@meta.data$barcodes <- rownames(macrophages@meta.data)
macrophages@meta.data <- macrophages@meta.data %>%
  left_join(
    md_meta %>% select(-project),  # remove duplicate column
    by = "sample"
  )
rownames(macrophages@meta.data) <- macrophages@meta.data$barcodes
DefaultAssay(macrophages) <- "RNA"

# ---- DE of interesting clusters ----
mac_c3 <- subset(macrophages,
                 subset = seurat_clusters == 3)
Idents(mac_c3) <- "aging"
deg_c3 <- FindMarkers(
  mac_c3,
  ident.1 = "aged",
  ident.2 = "young",
  test.use = "MAST",                       # or "wilcox"
  min.pct = 0.1,
  logfc.threshold = 0,                     # keep all; filter later (e.g., abs(log2FC)>=0.25)
  verbose = FALSE
)
# missing robust pseudobulk DE
