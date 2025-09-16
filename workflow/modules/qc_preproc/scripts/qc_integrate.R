library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(reticulate)

csv_path       <- "/datos/sensence/emilio/scall/config/proyects_metadata_samples_fastqs.csv"
root_obj       <- "/datos/sensence/emilio/scall/results/qc_preproc"
out_root       <- "/datos/sensence/emilio/scall/results/qc_preproc/integrated"

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

proj_col <- df$proyect
samp_col <- df$ident_sample

preintegrated_dir_for <- function(project, sample) {
  file.path(root_obj, project, sample, paste0(sample, "_20_qc.rds"))
}

objs <- list()
for (i in seq_len(nrow(df))) {
  project <- proj_col[i]; sample <- samp_col[i]
  path <- preintegrated_dir_for(project, sample)
  if (!file.exists(path)) { warning("Missing: ", path); next }
  obj <- readRDS(path)
  
  # add metadata (either assignment or AddMetaData; both are valid)
  obj$project <- project
  obj$sample  <- sample  # AddMetaData(obj, metadata = data.frame(project=project, sample=sample)) is also fine
  
  # ensure unique barcodes across samples
  obj <- RenameCells(obj, add.cell.id = sample)
  message("loading:",path)
  
  objs[[sample]] <- obj
}

message("Done loading")

# Merge them in a big seurat file
invivo.big <- merge(unlist(objs[1]), y = unlist(objs[-1]), project = "liver_hca")
# Separate each one in layers by sample (could change it to proyect)
invivo.big[["RNA"]] <- split(invivo.big[["RNA"]], f = invivo.big$sample)

# I dunno if we got to qc and preproc all over again

# Integrate using scVI since its the most robust method as of the current date 
obj <- IntegrateLayers(
  object = invivo.big, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "../miniconda3/envs/scvi-env", verbose = FALSE
)







