suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(patchwork)
})

# Figure creation for qc_preproc

in_vitro_seurat <- LoadSeuratRds("/datos/sensence/emilio/scall/results/qc_preproc/tecg_in_vitro.rds")
teiv <- LoadSeuratRds("/datos/sensence/emilio/scall/results/qc_preproc/te_in_vitro.rds")

in_vitro_seurat <- JoinLayers(in_vitro_seurat, assay = "RNA")

te <- rownames(teiv@assays[["RNA"]]) #TE list
var_genes <- VariableFeatures(in_vitro_seurat, nfeatures = 2000) #variable genes overall
var_te <- intersect(te, var_genes) #variable TE overall

umap_samp <- DimPlot(in_vitro_seurat, reduction = "harmony_umap", group.by = c("orig.ident"))
meta <- teiv@meta.data

# reshape to long format for ggplot
meta_long <- meta %>%
  select(orig.ident, nFeature_RNA, nCount_RNA) %>%
  pivot_longer(cols = c(nFeature_RNA, nCount_RNA),
               names_to = "Metric",
               values_to = "Value") %>%
  group_by(orig.ident) %>%
  mutate(median_val = median(Value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(orig.ident = reorder(orig.ident, median_val))


# violin plot
ggplot(meta_long, aes(x = orig.ident, y = Value, fill = orig.ident)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.3) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


stat_compare_means(
  method = "t.test",
  label = "p.signif",
  comparisons = list(
    c("control_DMSO", "etoposide_50m_10d"),
    c("control_DMSO", "irs_10gy")
  )
)

###############################

in_vitro_seurat$annotation=case_when(
  in_vitro_seurat$harmony_clusters == 7~"p16_high_ncount_low",
  
  in_vitro_seurat$harmony_clusters %in% c(2, 8 , 0, 1, 3)~"p16_high_ncount_high",
  
  in_vitro_seurat$harmony_clusters == 9~"p16_low_ncount_high",
  
  in_vitro_seurat$harmony_clusters %in% c(4, 6,5, 10)~"proliferating",
  
  TRUE~"others"
)


DimPlot(in_vitro_seurat, group.by ="annotation" )

in_vitro_seurat <- ScaleData(in_vitro_seurat, assay = "TE")
in_vitro_seurat <- RunPCA(in_vitro_seurat, npcs = 50, assay = "TE")
in_vitro_seurat = FindVariableFeatures(in_vitro_seurat, assay = "TE")
in_vitro_seurat[["TE"]] = split(in_vitro_seurat[["TE"]], in_vitro_seurat$orig.ident)
in_vitro_seurat = IntegrateLayers(in_vitro_seurat, HarmonyIntegration, new.reduction = "harmony", assay = "TE")
in_vitro_seurat = FindNeighbors(in_vitro_seurat, reduction = "harmony", dims = 1:30, assay = "TE")
in_vitro_seurat = FindClusters(in_vitro_seurat, resolution = 0.5, cluster.name = "harmony_clusters", assay = "TE")
in_vitro_seurat = RunUMAP(in_vitro_seurat, reduction = "harmony", dims = 1:30, reduction.name = "harmony_umap_TE", assay = "TE")
umap_samp <- DimPlot(in_vitro_seurat, reduction = "harmony_umap", group.by = c("orig.ident"))
umap_samp

teiv <- CreateSeuratObject(
  counts = in_vitro_TE@layers[["counts"]],   # raw counts
  assay = "RNA"
)
in_vitro_seurat = NormalizeData(teiv)
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

################################
#remotes::install_github('immunogenomics/presto')



all(colnames(in_vitro_seurat)==colnames(teiv))



in_vitro_seurat$annotation_TE=teiv$annotation_TE


DimPlot(in_vitro_seurat, group.by ="annotation_TE" )

FeaturePlot(in_vitro_seurat,"CCND1",min.cutoff = "q5",max.cutoff = "q95")

FeaturePlot(in_vitro_seurat,"NFIA",min.cutoff = "q5",max.cutoff = "q95")
FeaturePlot(in_vitro_seurat,"MKI67",min.cutoff = "q5",max.cutoff = "q95")
FeaturePlot(in_vitro_seurat,"NDUFA3",min.cutoff = "q5",max.cutoff = "q95")


in_vitro_seurat=JoinLayers(in_vitro_seurat)

markers_green_vs_red_RNA=FindMarkers(in_vitro_seurat, ident.1="cell_selector_right", ident.2="cell_selector_left", group.by ="annotation_TE" )

VariableFeatures(in_vitro_seurat)

markers_green_vs_red_RNA_filtered=markers_green_vs_red_RNA %>% filter(p_val_adj < 0.05, avg_log2FC> 0, pct.1>0.8)
