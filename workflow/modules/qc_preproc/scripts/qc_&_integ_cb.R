suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
})

#Variables
csv.path <- "/datos/sensence/emilio/liver_sc/fastq_prepros/config/proyects_metadata_samples_SC_fastqs.csv"
out.root <- "/datos/sensence/emilio/scall/results/qc_preproc"
combined.root <- "/datos/sensence/emilio/scall/results/te"

############# Load objects ################
meta <- read.csv(csv.path, stringsAsFactors = FALSE)
meta_unique <- unique(meta[, c("proyect", "ident_sample")])

seurat_objects <- list()
for (i in seq_len(nrow(meta_unique))) {
  project <- meta_unique$proyect[i]
  sample  <- meta_unique$ident_sample[i]
  
  sample_path <- file.path(combined.root, project, sample, "combined_corrected")
  message("Loading sample: ", project, " / ", sample)
  
  # Read TE and GENE layers
  combined_layers <- Read10X(sample_path)
  combined_counts <- rbind(combined_layers[["Gene"]], combined_layers[["TE"]])
  
  # Make Seurat object
  sample_seurat <- CreateSeuratObject(
    counts = combined_counts,
    project = paste(project, sample, sep = "_"),
    min.cells = 3,
    min.features = 50
  )
  
  # Convert to Seurat v5 assay
  sample_seurat[["RNA"]] <- as(sample_seurat[["RNA"]], "Assay5")
  # Add the Gene expresion and TE layers from the combined_layers obj
  sample_seurat[["CG"]] <- CreateAssayObject(counts = combined_layers[["Gene"]])
  sample_seurat[["CG"]] <- as(sample_seurat[["CG"]], "Assay5")
  sample_seurat[["TE"]] <- CreateAssayObject(counts = combined_layers[["TE"]])
  sample_seurat[["TE"]] <- as(sample_seurat[["TE"]], "Assay5")
  
  # Store in list
  seurat_objects[[paste(project, sample, sep = "_")]] <- sample_seurat
  
}

in_vitro_seurat = Reduce(f = merge, seurat_objects)

############# QC ############

in_vitro_seurat = NormalizeData(in_vitro_seurat)
# join the layers to calculate mitochondrial percentage 
in_vitro_seurat = JoinLayers(in_vitro_seurat)
in_vitro_seurat = PercentageFeatureSet(in_vitro_seurat, pattern = "^MT-", col.name = "percent.mt")
in_vitro_seurat = subset(in_vitro_seurat, percent.mt < 30)
# split the layers again. this is required for harmony
in_vitro_seurat[["RNA"]] = split(in_vitro_seurat[["RNA"]], in_vitro_seurat$orig.ident)
in_vitro_seurat = FindVariableFeatures(in_vitro_seurat)
in_vitro_seurat = ScaleData(in_vitro_seurat)
in_vitro_seurat = RunPCA(in_vitro_seurat, npcs = 30)
in_vitro_seurat = IntegrateLayers(in_vitro_seurat, HarmonyIntegration, new.reduction = "harmony")
in_vitro_seurat = FindNeighbors(in_vitro_seurat, reduction = "harmony", dims = 1:30)
in_vitro_seurat = FindClusters(in_vitro_seurat, resolution = 0.5, cluster.name = "harmony_clusters")
in_vitro_seurat = RunUMAP(in_vitro_seurat, reduction = "harmony", dims = 1:30, reduction.name = "harmony_umap")





