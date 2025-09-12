suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(R.utils)
  library(dplyr)
  library(patchwork)
})

# ---- USER PATHS ----
csv_path       <- "/datos/sensence/emilio/scall/config/proyects_metadata_samples_fastqs.csv"
root_counts    <- "/datos/sensence/emilio/liver_sc/fastq_prepros/output/counts"
out_root       <- "/datos/sensence/emilio/scall/results/qc_preproc"

# ---- QC PARAMS ----
layer          <- "filtered"   # or "raw"
mito_regex     <- "^MT-"       # Human; use "^mt-" for mouse
min_genes      <- 200
max_mt_pct     <- 14
convert_h5ad   <- TRUE
i              <- 1

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# ---- READ CSV (auto-detect project/sample columns) ----
df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

proj_col <- df$proyect
samp_col <- df$ident_sample

# ---- HELPERS ----
mtx_dir_for <- function(project, sample, layer) {
  file.path(root_counts, project, paste0(sample, "_matrix_Solo.out"), "Gene", layer)
}

project <- proj_col[i]
sample  <- samp_col[i]
mtx_dir <- mtx_dir_for(project, sample, layer)

if (!dir.exists(mtx_dir)) {
  message(sprintf("SKIP: %s/%s -- missing %s", project, sample, mtx_dir))
  next
}
message(sprintf("[QC] %s / %s layer=%s", project, sample, layer))

#see if it is gzipped as read10x needs it 

is_gz <- function(path) {
  if (!file.exists(path)) stop("No such file: ", path)
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  sig <- readBin(con, "raw", 2L)
  length(sig) == 2L && as.integer(sig[1]) == 0x1f && as.integer(sig[2]) == 0x8b
}

# Compress in-place with system gzip if needed; returns gz path
gzip_if_needed <- function(path) {
  stopifnot(file.exists(path), !dir.exists(path))
  if (is_gz(path)) {
    if (!grepl("\\.gz$", path)) {  # gz content but no .gz suffix → rename
      dest <- paste0(path, ".gz"); file.rename(path, dest); return(dest)
    }
    return(path)                   # already gz and named *.gz
  }
  # plain file; if name ends with .gz but content isn't gz, fix name first
  if (grepl("\\.gz$", path)) { base <- sub("\\.gz$", "", path); file.rename(path, base); path <- base }
  system2("gzip", c("-f", shQuote(path)), stdout = NULL, stderr = NULL)  # -f = force overwrite if *.gz exists
  paste0(path, ".gz")
}

gz_paths <- sapply(
  c("matrix.mtx","barcodes.tsv","features.tsv"),
  function(f) { f <- file.path(mtx_dir, f); if (file.exists(f)) gzip_if_needed(f) },
  USE.NAMES = TRUE
)

#Load it to seurat since regex wont work now with Read10x function
pick <- function(d, a, b) if (file.exists(file.path(d,a))) file.path(d,a) else file.path(d,b)
mtx  <- pick(mtx_dir, "matrix.mtx.gz",  "matrix.mtx")
bc   <- pick(mtx_dir, "barcodes.tsv.gz","barcodes.tsv")
feat <- {
  f <- pick(mtx_dir, "features.tsv.gz","features.tsv")
  if (!file.exists(f)) f <- pick(mtx_dir, "genes.tsv.gz","genes.tsv")
  f
}

# Count columns in the first line of features to choose the right name column
feat_cols <- length(strsplit(readLines(feat, n = 1), "\t", fixed = TRUE)[[1]])
gene_col  <- if (feat_cols >= 2) 2 else 1

#load seurat obj
counts <- ReadMtx(mtx = mtx, cells = bc, features = feat, feature.column = gene_col)
obj <- CreateSeuratObject(counts)

#Start trimming the outliers
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA >= min_genes & percent.mt < max_mt_pct)

#normalize data
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature selection
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

#Scale the data (here we can also regress mit percentage or cell cycle)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)

# Linear Dimensional reduction
obj <- RunPCA(obj, features = VariableFeatures(object = obj))

#Choose how many dimensionalities we want
#on process since i need to update all the f ing biocon 

#clustering 
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)

# Non linear Dimensional reduction
obj <- RunUMAP(obj, dims = 1:10)

proj_out <- file.path(out_root, project, sample)
dir.create(proj_out, showWarnings = FALSE, recursive = TRUE)

h5s_path <- file.path(proj_out, paste0(sample, "_20_qc.h5seurat"))
SaveH5Seurat(obj, filename = h5s_path, overwrite = TRUE)

if (isTRUE(convert_h5ad)) {
  SeuratDisk::Convert(h5s_path, dest = "h5ad", overwrite = TRUE)