te_seurat <- in_vitro_seurat
Assays(te_seurat)

DefaultAssay(te_seurat) = "TE"
DefaultAssay(te_seurat)

te_seurat = NormalizeData(te_seurat)

# split the layers again. this is required for harmony
te_seurat[["TE"]] = split(te_seurat[["TE"]], te_seurat$orig.ident)

te_seurat = FindVariableFeatures(te_seurat)
te_seurat = ScaleData(te_seurat)
te_seurat = RunPCA(te_seurat, npcs = 30)

te_seurat = IntegrateLayers(te_seurat, HarmonyIntegration, new.reduction = "harmony")
te_seurat = FindNeighbors(te_seurat, reduction = "harmony", dims = 1:30)
te_seurat = RunUMAP(
  te_seurat,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "harmony_umap"
)



te_seurat = FindClusters(te_seurat, resolution = 0.5, cluster.name = "harmony_clusters")

DimPlot(te_seurat)

cell_selector_left = CellSelector(DimPlot(te_seurat, group.by = "annotation"), te_seurat)
cell_selector_left = Idents(cell_selector_left)
cell_selector_left = names(cell_selector_left)[cell_selector_left == "SelectedCells"]


cell_selector_right = CellSelector(DimPlot(te_seurat, group.by = "annotation"), te_seurat)
cell_selector_right = Idents(cell_selector_right)
cell_selector_right = names(cell_selector_right)[cell_selector_right == "SelectedCells"]

length(cell_selector_right)
length(cell_selector_left)

###############################
te_seurat$annotation_TE = case_when(
  colnames(te_seurat) %in% cell_selector_left ~ "cell_selector_left",
  colnames(te_seurat) %in% cell_selector_right ~ "cell_selector_right",
  TRUE ~ "others"
)



DimPlot(te_seurat, group.by = "annotation")


te_seurat@meta.data  %>%  group_by(annotation_TE)  %>%  mutate(annotation_TE_total =
                                                                 n()) %>% ungroup() %>% group_by(annotation_TE, annotation)  %>%
  summarize(percentage = n() * 100 / annotation_TE_total) %>%  distinct() %>% pivot_wider(names_from = annotation, values_from =
                                                                                            percentage)



te_seurat@meta.data  %>%  group_by(annotation_TE)  %>%  mutate(annotation_TE_total =
                                                                 n()) %>% ungroup() %>% group_by(annotation_TE, annotation)  %>%
  summarize(percentage = n() * 100 / annotation_TE_total) %>%  distinct()  %>%
  ggplot(aes(x = annotation_TE, y = percentage, fill = annotation)) + geom_col()


###############################
te_seurat = JoinLayers(te_seurat)
markers_stressed_vs_controls_TE = FindMarkers(te_seurat,
                                              group.by = "annotation_TE",
                                              ident.1 = "cell_selector_right",
                                              ident.2 = "cell_selector_left")



markers_stressed_vs_controls_TE_family = markers_stressed_vs_controls_TE %>%
  filter(p_val_adj < 0.05) %>%
  mutate(family = stringr::str_extract(rownames(.), pattern = "^[A-Z]\\w+")) %>%
  mutate(inverted_p_value = -log10(ifelse(
    p_val_adj == 0, .Machine$double.xmin , p_val_adj
  ))) %>%
  mutate(weighted_log2fc = avg_log2FC * inverted_p_value) %>%
  ##########
group_by(family) %>%
  summarize(family_fold_change = mean(weighted_log2fc))



markers_stressed_vs_controls_TE_family %>%
  #
  ggplot(aes(x = family_fold_change, y = reorder(family, family_fold_change))) +
  geom_col()+
    labs(x="Pooled weighted fold change per family",y="TE Family") + theme_light()




###############################
te_seurat$orig.ident %>% table

DimPlot(te_seurat,group.by =  "orig.ident")
