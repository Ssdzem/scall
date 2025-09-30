suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(patchwork)
})

# Variables
csv_path <- "/datos/sensence/emilio/liver_sc/fastq_prepros/config/proyects_metadata_samples_SC_fastqs.csv"
out_root <- "/datos/sensence/emilio/scall/results/qc_preproc"
combined_root <- "/datos/sensence/emilio/scall/results/te"

############# Load objects ################
meta <- read.csv(csv_path, stringsAsFactors = FALSE)
meta_unique <- unique(meta[, c("proyect", "ident_sample")])

# for the combined Canonical genes (CG) + transposable elements (TE)
seurat_objects_tecg <- list()
# only TE
seurat_objects_te <- list()
# only CG
seurat_objects_cg <- list()

for (i in seq_len(nrow(meta_unique))) {
  project <- meta_unique$proyect[i]
  sample <- meta_unique$ident_sample[i]
  sample_path <- file.path(combined_root, project, sample, "combined_corrected")
  message("Loading sample: ", project, " / ", sample)
  # Read TE and CG layers
  combined_layers <- Read10X(sample_path)
  combined_counts <- rbind(combined_layers[["Gene"]], combined_layers[["TE"]])
  # Make Seurat object
  sample_seurat_tecg <- CreateSeuratObject(
    counts = combined_counts,
    project = paste(project, sample, sep = "_")
  )
  sample_seurat_te <- CreateSeuratObject(
    counts = combined_layers[["TE"]],
    project = paste(project, sample, sep = "_")
  )
  sample_seurat_cg <- CreateSeuratObject(
    counts = combined_layers[["Gene"]],
    project = paste(project, sample, sep = "_")
  )
  # Convert to Seurat v5 assay
  sample_seurat_tecg[["RNA"]] <- as(sample_seurat_tecg[["RNA"]], "Assay5")
  sample_seurat_te[["RNA"]] <- as(sample_seurat_te[["RNA"]], "Assay5")
  sample_seurat_cg[["RNA"]] <- as(sample_seurat_cg[["RNA"]], "Assay5")
  # Store in its respective list
  seurat_objects_tecg[[paste(project, sample, sep = "_")]] <- sample_seurat_tecg
  seurat_objects_te[[paste(project, sample, sep = "_")]] <- sample_seurat_te
  seurat_objects_cg[[paste(project, sample, sep = "_")]] <- sample_seurat_cg
}

in_vitro_seurat_tecg <- Reduce(f = merge, seurat_objects_tecg)
in_vitro_seurat_te <- Reduce(f = merge, seurat_objects_te)
in_vitro_seurat_cg <- Reduce(f = merge, seurat_objects_cg)

############# QC ############
qc_pipeline <- function(obj) {
  obj <- NormalizeData(obj)
  obj <- JoinLayers(obj, assay = "RNA")
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
  obj <- subset(obj, subset = percent.mt < 30)
  obj[["RNA"]] <- split(obj[["RNA"]], obj$orig.ident)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 50)
  obj <- IntegrateLayers(obj, HarmonyIntegration, new.reduction = "harmony")
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5, cluster.name = "harmony_clusters")
  obj <- RunUMAP(obj, reduction = "harmony",
                 dims = 1:50, reduction.name = "harmony_umap")
  obj
}

objs_in <- list(
  tecg = in_vitro_seurat_tecg,
  te = in_vitro_seurat_te,
  cg = in_vitro_seurat_cg
)

objs_out <- lapply(objs_in, qc_pipeline)

in_vitro_seurat_tecg <- objs_out$tecg
in_vitro_seurat_te <- objs_out$te
in_vitro_seurat_cg <- objs_out$cg

saveRDS(in_vitro_seurat_tecg, file = file.path(out_root, "tecg_in_vitro.rds"))
saveRDS(in_vitro_seurat_te, file = file.path(out_root, "te_in_vitro.rds"))
saveRDS(in_vitro_seurat_cg, file = file.path(out_root, "cg_in_vitro.rds"))