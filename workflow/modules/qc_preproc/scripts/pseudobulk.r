library(Seurat)
library(tidyverse)

tecg_path <- "/datos/sensence/emilio/scall/results/qc_preproc/tecg_in_vitro.rds"
tecg <- LoadSeuratRds(tecg_path)

pseudo_tecg <- AggregateExpression(tecg,
                                   assays = "RNA",
                                   return.seurat = TRUE,
                                   group.by = c("orig.ident"))
Cells(pseudo_tecg)
Idents(pseudo_tecg)

controls <- c("control-DMSO", "control-PDL24")
pseudo_tecg$cellstate <- ifelse(Idents(pseudo_tecg) %in% controls,
                                "proliferative", "senescent")
pseudo_tecg$cellstate <- factor(pseudo_tecg$cellstate,
                                levels = c("proliferative","senescent"))
table(Idents(pseudo_tecg), pseudo_tecg$cellstate)
Idents(pseudo_tecg) <- "cellstate"

te_path <- "/datos/sensence/emilio/scall/results/qc_preproc/te_in_vitro.rds"
te <- LoadSeuratRds(te_path)

te <- rownames(te)

DE_pseudo <- FindMarkers(object = pseudo_tecg,
                         ident.1 = "proliferative",
                         ident.2 = "senescent",
                         min.cells.group = 2,
                         min.pct = 0,
                         test.use = "DESeq2")

DE_te_tbl <- DE_pseudo %>%
  rownames_to_column(var = "TECG") %>%
  filter(TECG %in% te) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(avg_log2FC) %>%
  mutate(
    direction = ifelse(avg_log2FC > 0, "Up in proliferative", "Up in senescent")
  )

DE_te_tbl$TECG <- factor(DE_te_tbl$TECG, levels = DE_te_tbl$TECG)

ggplot(DE_te_tbl, aes(x = avg_log2FC, y = TECG, fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Up in proliferative" = "steelblue",
                               "Up in senescent" = "firebrick")) +
  labs(
    x = "log2 Fold Change (proliferative / senescent)",
    y = "TE (FDR < 0.05)",
    title = "Differentially Expressed TEs in Pseudobulk in vitro dataset"
  ) +
  theme_bw()
