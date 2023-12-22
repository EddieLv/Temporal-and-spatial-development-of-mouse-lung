#####################################
# seed may cause slightly different!
#####################################

set.seed(666)
source("/media/biogenger/D/scripts/Spatial-Transcriptome/R/initialize.genger_env.R")
source("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/up.markers.R")
source("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/down.markers.R")
initialize.genger_env()

### 0.Initialize ###
# Load meta.info
meta.info = readxl::read_xlsx("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/meta.info.xlsx") %>%
  filter(time %in% c("E12.5", "E15.5", "E16.5", "E17.5", "E18.5", "P0", "P3") & !(name %in% c("E16.5_slice18", "E18.5_slice31")))
srat.merge = readRDS("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/1-select_feature/cache-select_parameter/srat.merge.raw.pre.rds")
srat.merge = RunUMAP(srat.merge, reduction = "pca", dims = 1:36, umap.method = "uwot", metric = "cosine", min.dist = 0.01, n.neighbors = 20, seed.use = 666, return.model = T)
################################################################################################################################
library(scCustomize)
### 1.Dimensionality-reduction & Clustering ###
res.list = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.3, 1.4, 1.5)
library(clustree)
clustree(srat.merge, prefix = "SCT_snn_res.", node_colour = "grey70", show_axis = T)
wrap_plots(map(res.list, ~DimPlot_scCustom(srat.merge, group.by = paste0("SCT_snn_res.", .), label = T, repel = F)), ncol = 4) & theme(aspect.ratio = 1)

FeaturePlot_scCustom(srat.merge, features = c("Scgb3a2", "Scgb1a1", "Foxj1", "Deup1", "Tgfbi", "Wnt5a", "Ager"), colors_use = viridis(n = 100, option = "D"), order = T) &
  theme(aspect.ratio = 1) & labs(subtitle = "sct_normalized_data")
DefaultAssay(srat.merge) = "MAGIC_SCT"
FeaturePlot_scCustom(srat.merge, features = c("Scgb3a2", "Scgb1a1", "Foxj1", "Deup1", "Tgfbi", "Wnt5a", "Ager"), colors_use = viridis(n = 100, option = "D"), order = T) &
  theme(aspect.ratio = 1) & labs(subtitle = "denoised")
DefaultAssay(srat.merge) = "SCT"

# [PCs] and [Resolution] are determined!
res_use = "SCT_snn_res.1.5"
srat.merge[["seurat_clusters"]] = srat.merge[[res_use]]
srat.merge = SetIdent(srat.merge, value = "seurat_clusters")
p1 = DimPlot_scCustom(srat.merge, group.by = "orig.ident", figure_plot = T)
p2 = DimPlot_scCustom(srat.merge, group.by = "ind_anno_ctp3", label = T, repel = T, colors_use = pal, figure_plot = T)
p3 = DimPlot_scCustom(srat.merge, group.by = "seurat_clusters", label = T, repel = T, figure_plot = T)
wrap_plots(p1, p2, p3, ncol = 3) & theme(aspect.ratio = 1)

# Immune
Stacked_VlnPlot(srat.merge, features = gene.sets.up$immune$immune, x_lab_rotate = T)
# Epithelium
Stacked_VlnPlot(srat.merge, features = gene.sets.up$epith$epith, x_lab_rotate = T)
Stacked_VlnPlot(srat.merge, features = c("Foxj1", "Scgb1a1", "Cdk8", "Il31ra", "Trim30a", "Ager"), x_lab_rotate = T)
# Endothelium
Stacked_VlnPlot(srat.merge, features = gene.sets.up$endo$endo, x_lab_rotate = T)
# Mesenchymal
Stacked_VlnPlot(srat.merge, features = gene.sets.up$mesen$mesen, x_lab_rotate = T)
Stacked_VlnPlot(srat.merge, features = c("Tgfbi", "Wnt5a", "Sftpb", "Cxcl15"), x_lab_rotate = T)
# Unknown
Stacked_VlnPlot(srat.merge, features = c("Cdk8", "Il31ra", "Trim30a", "Trim30b", "Trim30c", "Trim12c"), x_lab_rotate = T)

cluster.heat(srat.merge, assay = "SCT", slot = "data", cluster_col = "celltype.merge")
cluster.tree(srat.merge, assay = "SCT", reduction = "pca", dims = 1:36)

Stacked_VlnPlot(srat.merge, features = c("Samd5", "Vwf", "Vegfc", "Jam2", "Ephb4"), x_lab_rotate = T) # Venous

### 2.Refine Cluster ###
srat.merge$split_ctp = as.character(srat.merge$seurat_clusters)
DimPlot_scCustom(srat.merge, group.by = "seurat_clusters")

refine.samples = c("Mesenchymal", "Alveolar", "Bronchi", "Endothelium", "Unknown")
for (sample in refine.samples) {
  gc()
  message(paste("Start", sample))
  srat.sub = readRDS(paste0("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/2-split_annotation/cache-split_annotation/srat.sub.", sample, ".anno.rds"))
  srat.merge$split_ctp[Cells(srat.sub)] = srat.sub$small.type
}
DimPlot_scCustom(srat.merge, group.by = "split_ctp", label = T, repel = T)
####### Rename ######
srat.merge$split_ctp[srat.merge$split_ctp == "AT12"] = "Alveolar.Mixture"
srat.merge$split_ctp[srat.merge$split_ctp == "Ciliated"] = "Ciliated.like"
srat.merge$split_ctp[srat.merge$split_ctp == "Secretory"] = "Secretory.like"
srat.merge$split_ctp[srat.merge$split_ctp == "PNEC"] = "PNEC.like"
###########Complete Name ###########
srat.merge$full.name = as.character(srat.merge$split_ctp)
srat.merge$full.name[srat.merge$full.name == "Sox2.pro"] = "Sox2+ epithelial progenitor"
srat.merge$full.name[srat.merge$full.name == "Ciliated.like"] = "Ciliated cell"
srat.merge$full.name[srat.merge$full.name == "Secretory.like"] = "Secretory cell"
srat.merge$full.name[srat.merge$full.name == "PNEC.like"] = "Pulmonary neuroendocrine cell (PNEC)"
srat.merge$full.name[srat.merge$full.name == "Sox9.pro"] = "Sox9+ epithelial progenitor"
srat.merge$full.name[srat.merge$full.name == "AT1.pre"] = "AT1 precursor"
srat.merge$full.name[srat.merge$full.name == "AT2.pre"] = "AT2 precursor"
srat.merge$full.name[srat.merge$full.name == "Alveolar.Mixture"] = "Mature alveolar cell"
srat.merge$full.name[srat.merge$full.name == "Matrix.FB"] = "Matrix fibroblast"
srat.merge$full.name[srat.merge$full.name == "PMP"] = "Proliferating matrix fibroblast"
srat.merge$full.name[srat.merge$full.name == "myo.FB"] = "Myofibroblast"
srat.merge$full.name[srat.merge$full.name == "ASMC.pre"] = "ASMC precursor"
srat.merge$full.name[srat.merge$full.name == "Il31ra.FB"] = "Il31ra+ mesenchymal cells"
srat.merge$full.name[srat.merge$full.name == "EC.pro"] = "Endothelium progenitor"
srat.merge$full.name[srat.merge$full.name == "AEC"] = "Arterial endothelium"
srat.merge$full.name[srat.merge$full.name == "VEC"] = "Venous endothelium"
srat.merge$full.name[srat.merge$full.name == "LEC"] = "Lymphatic endothelium"
srat.merge$full.name[srat.merge$full.name == "VSMC.pre"] = "VSMC precursor"
srat.merge$full.name[srat.merge$full.name == "ASMC"] = "Airway smooth muscle cell (ASMC)"
srat.merge$full.name[srat.merge$full.name == "VSMC"] = "Vascular smooth muscle cell (VSMC)"
srat.merge$full.name[srat.merge$full.name == "adventitial.FB"] = "Adventitial fibroblast"
srat.merge$full.name = factor(srat.merge$full.name, levels = c("Sox2+ epithelial progenitor", "Ciliated cell", "Secretory cell", "Pulmonary neuroendocrine cell (PNEC)",
                                                               "Sox9+ epithelial progenitor", "AT1 precursor", "AT2 precursor", "Mature alveolar cell",
                                                               "Proliferating matrix fibroblast", "Matrix fibroblast", "Myofibroblast", "ASMC precursor", "Airway smooth muscle cell (ASMC)", "Mesothelium", "Il31ra+ mesenchymal cells",
                                                               "Endothelium progenitor", "Arterial endothelium", "Venous endothelium", "Lymphatic endothelium",
                                                               "VSMC precursor", "Vascular smooth muscle cell (VSMC)", "Adventitial fibroblast", "Myocyte",
                                                               "Chondrocyte"))
srat.merge$full.name.number = as.character(srat.merge$full.name)
for (i in 1:length(levels(srat.merge$full.name))) {
  srat.merge$full.name.number[srat.merge$full.name == levels(srat.merge$full.name)[i]] = paste0(i, " ", levels(srat.merge$full.name)[i])
}
srat.merge$full.name.number = factor(as.character(srat.merge$full.name.number), levels = str_sort(unique(srat.merge$full.name.number), numeric = T))

DimPlot_scCustom(srat.merge, group.by = "full.name")
# srat.merge$full.name.comb = "oops"
# for (i in 1:length(levels(srat.merge$full.name))) {
#   srat.merge$full.name.comb[srat.merge$full.name == levels(srat.merge$full.name)[i]] = paste0(i, ".", levels(srat.merge$full.name)[i])
# }
# srat.merge$full.name.comb = factor(as.character(srat.merge$full.name.comb), levels = str_sort(unique(srat.merge$full.name.comb), numeric = T))
# DimPlot_scCustom(srat.merge, group.by = "full.name.comb")
#####################
srat.merge = SetIdent(srat.merge, value = "split_ctp")
DimPlot_scCustom(srat.merge, group.by = "split_ctp", label = T, repel = T)

library(MetaNeighbor)
srat.merge.sce = as.SingleCellExperiment(srat.merge)
var.genes = variableGenes(dat = srat.merge.sce, exp_labels = srat.merge$orig.ident)
head(var.genes)
aurocs = MetaNeighborUS(var_genes = var.genes, dat = srat.merge.sce, study_id = srat.merge$orig.ident, cell_type = srat.merge$split_ctp, fast_version = T)
MetaNeighbor::plotHeatmap(aurocs, cex = 1)
rm(srat.merge.sce); gc()

pdf("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/split_ctp.similarity.pdf", width = 8*1.1, height = 8)
cluster.tree(srat.merge, assay = "SCT", reduction = "pca", dims = 1:36)
cluster.heat(srat.merge, assay = "SCT", slot = "data", cluster_col = "split_ctp")
dev.off()

srat.merge$merge_ctp2 = "good"
srat.merge$merge_ctp2[WhichCells(srat.merge, idents = c("Matrix.FB", "Il31ra.FB", "Mesothelium", "myo.FB", "PMP", "ASMC", "ASMC.pre"))] = "Mesenchyme"
srat.merge$merge_ctp2[WhichCells(srat.merge, idents = c("AT1.pre", "AT2.pre", "Alveolar.Mixture", "Sox9.pro"))] = "Alveolar Epithlium"
srat.merge$merge_ctp2[WhichCells(srat.merge, idents = c("Sox2.pro", "Ciliated.like", "Secretory.like", "PNEC.like"))] = "Bronchi Epithlium"
srat.merge$merge_ctp2[WhichCells(srat.merge, idents = c("AEC", "VEC", "LEC", "EC.pro"))] = "Vascular Endothelium"
srat.merge$merge_ctp2[WhichCells(srat.merge, idents = c("adventitial.FB", "VSMC.pre", "VSMC", "Myocyte"))] = "Vascular Fibroblast"
srat.merge$merge_ctp2[WhichCells(srat.merge, idents = c("Chondrocyte"))] = "Chondrocyte"
srat.merge = SetIdent(srat.merge, value = "merge_ctp2")
DimPlot_scCustom(srat.merge, group.by = "merge_ctp2", label = T, repel = T)
pdf("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/merge_ctp2.similarity.pdf", width = 8*1.1, height = 8)
cluster.tree(srat.merge, assay = "SCT", reduction = "pca", dims = 1:36)
cluster.heat(srat.merge, assay = "SCT", slot = "data", cluster_col = "merge_ctp2")
dev.off()

srat.merge$merge_ctp1 = "job"
srat.merge$merge_ctp1[WhichCells(srat.merge, idents = c("Mesenchyme"))] = "Mesenchyme"
srat.merge$merge_ctp1[WhichCells(srat.merge, idents = c("Alveolar Epithlium", "Bronchi Epithlium"))] = "Airway"
srat.merge$merge_ctp1[WhichCells(srat.merge, idents = c("Vascular Endothelium", "Vascular Fibroblast"))] = "Vascular"
srat.merge$merge_ctp1[WhichCells(srat.merge, idents = c("Chondrocyte"))] = "Chondrocyte"
srat.merge = SetIdent(srat.merge, value = "merge_ctp1")
DimPlot_scCustom(srat.merge, group.by = "merge_ctp1", label = T, repel = T)
pdf("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/merge_ctp1.similarity.pdf", width = 8*1.1, height = 8)
cluster.tree(srat.merge, assay = "SCT", reduction = "pca", dims = 1:36)
cluster.heat(srat.merge, assay = "SCT", slot = "data", cluster_col = "merge_ctp1")
dev.off()

srat.merge$split_ctp = factor(as.character(srat.merge$split_ctp), levels = c("Sox2.pro", "Ciliated.like", "Secretory.like", "PNEC.like",
                                                                             "Sox9.pro", "AT1.pre", "AT2.pre", "Alveolar.Mixture",
                                                                             "PMP", "Matrix.FB", "myo.FB", "ASMC.pre", "ASMC", "Mesothelium", "Il31ra.FB",
                                                                             "EC.pro", "AEC", "VEC", "LEC",
                                                                             "VSMC.pre", "VSMC", "adventitial.FB", "Myocyte",
                                                                             "Chondrocyte"))
srat.merge$merge_ctp2 = factor(as.character(srat.merge$merge_ctp2), levels = c("Alveolar Epithlium", "Bronchi Epithlium", "Mesenchyme", "Vascular Fibroblast", "Vascular Endothelium", "Chondrocyte"))
srat.merge$merge_ctp1 = factor(as.character(srat.merge$merge_ctp1), levels = c("Airway", "Mesenchyme", "Vascular", "Chondrocyte"))
srat.merge$orig.ident = factor(as.character(srat.merge$orig.ident), levels = str_sort(unique(srat.merge$orig.ident), numeric = T))
srat.merge$timepoint = factor(as.character(srat.merge$timepoint), levels = c("E12.5", "E15.5", "E16.5", "E17.5", "E18.5", "P0", "P3"))
# srat.merge$num.type = 666
# for (i in 1:length(levels(srat.merge$split_ctp))) {
#   srat.merge$num.type[srat.merge$split_ctp == levels(srat.merge$split_ctp)[i]] = i
# }
# srat.merge$num.type = factor(as.character(srat.merge$num.type), levels = str_sort(unique(srat.merge$num.type), numeric = T))
# srat.merge$comb.type = "oops"
# for (i in 1:length(levels(srat.merge$split_ctp))) {
#   srat.merge$comb.type[srat.merge$split_ctp == levels(srat.merge$split_ctp)[i]] = paste0(i, ".", levels(srat.merge$split_ctp)[i])
# }
# srat.merge$comb.type = factor(as.character(srat.merge$comb.type), levels = str_sort(unique(srat.merge$comb.type), numeric = T))

library(SCP)
p1 = CellDimPlot(srat.merge, group.by = "seurat_clusters", reduction = "UMAP", theme_use = "theme_blank", palcolor = pal[6:(6+length(unique(srat.merge$seurat_clusters)))],
                 show_stat = T, label = T, label_insitu = T, label_repel = T, label.size = 6, label_repulsion = 20, label_segment_color = "red", label_point_color = "black",
                 legend.position = "top", legend.direction = "horizontal") & theme(aspect.ratio = 1)
p2 = CellDimPlot(srat.merge, group.by = "split_ctp", reduction = "UMAP", theme_use = "theme_blank", palcolor = colors.lung.merge[levels(srat.merge$split_ctp)],
                 show_stat = T, label = T, label_insitu = T, label_repel = T, label.size = 6, label_repulsion = 20, label_segment_color = "red", label_point_color = "black",
                 legend.position = "top", legend.direction = "horizontal") & theme(aspect.ratio = 1)
p3 = CellDimPlot(srat.merge, group.by = "full.name", reduction = "UMAP", theme_use = "theme_blank", palcolor = colors.lung.merge[levels(srat.merge$split_ctp)],
                 show_stat = T, label = T, label_insitu = T, label_repel = T, label.size = 5, label_repulsion = 20, label.bg = NA, label.fg = NA, label_segment_color = "red", label_point_color = "black",
                 legend.position = "top", legend.direction = "horizontal") & theme(aspect.ratio = 1)
p4 = CellDimPlot(srat.merge, group.by = "full.name.number", reduction = "UMAP", theme_use = "theme_blank", palcolor = colors.lung.merge[levels(srat.merge$split_ctp)],
                 show_stat = T, label = T, label_insitu = T, label_repel = T, label.size = 5, label_repulsion = 20, label_segment_color = "red", label_point_color = "black",
                 legend.position = "top", legend.direction = "horizontal") & theme(aspect.ratio = 1)
srat.merge$tmp = factor(str_extract(srat.merge$full.name.number, "\\d+"), levels = 1:length(unique(srat.merge$full.name.number)))
p5 = CellDimPlot(srat.merge, group.by = "tmp", reduction = "UMAP", theme_use = "theme_blank", palcolor = colors.lung.merge[levels(srat.merge$split_ctp)],
                 show_stat = T, label = T, label_insitu = T, label_repel = T, label.size = 8, label_repulsion = 20, label_segment_color = "red", label_point_color = "black",
                 legend.position = "top", legend.direction = "horizontal") & theme(aspect.ratio = 1)
p6 = CellDimPlot(srat.merge, group.by = "merge_ctp2", reduction = "UMAP", theme_use = "theme_blank", palcolor = colors.mid.type[levels(srat.merge$merge_ctp2)],
                 show_stat = T, label = T, label_insitu = T, label_repel = T, label.size = 6, label_repulsion = 20, label_segment_color = "red", label_point_color = "black",
                 legend.position = "top", legend.direction = "horizontal") & theme(aspect.ratio = 1)
p7 = CellDimPlot(srat.merge, group.by = "merge_ctp1", reduction = "UMAP", theme_use = "theme_blank", palcolor = colors.big.type[levels(srat.merge$merge_ctp1)],
                 show_stat = T, label = T, label_insitu = T, label_repel = T, label.size = 6, label_repulsion = 20, label_segment_color = "red", label_point_color = "black",
                 legend.position = "top", legend.direction = "horizontal") & theme(aspect.ratio = 1)
p8 = CellDimPlot(srat.merge, group.by = "orig.ident", reduction = "UMAP", theme_use = "theme_blank", palcolor = pal[6:(6+length(unique(srat.merge$orig.ident)))],
                 show_stat = T, label = T, label_insitu = T, label_repel = T, label.size = 6, label_repulsion = 20, label_segment_color = "red", label_point_color = "black",
                 legend.position = "top", legend.direction = "horizontal") & theme(aspect.ratio = 1)
p9 = CellDimPlot(srat.merge, group.by = "timepoint", reduction = "UMAP", theme_use = "theme_blank", palcolor = colors.timepoint[levels(srat.merge$timepoint)],
                 show_stat = T, label = F, label_insitu = T, label_repel = T, label.size = 6, label_repulsion = 20, label_segment_color = "red", label_point_color = "black",
                 legend.position = "top", legend.direction = "horizontal") & theme(aspect.ratio = 1)
pdf("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/annotation.umap.pdf", width = 8*1.7, height = 8*1.7)
p1
p2
p3
p4
p5
p6
p7
p8
p9
dev.off()

plots = CellDimPlot(srat.merge, group.by = "full.name", split.by = "full.name", reduction = "UMAP", palcolor = colors.lung.merge[levels(srat.merge$split_ctp)],
                    label = T, label_insitu = T, label_repel = T, label.size = 6, label_repulsion = 80, label_segment_color = "red", label_point_color = "black",
                    show_stat = T, cells.highlight = T, theme_use = "theme_blank", legend.position = "none", combine = F)
pdf("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/each.full.name.umap.pdf", width = 8, height = 8)
walk(plots, ~print(.))
dev.off()

plots = CellDimPlot(srat.merge, group.by = "timepoint", split.by = "timepoint", reduction = "UMAP", palcolor = colors.timepoint[levels(srat.merge$timepoint)],
                    label = F, label_insitu = T, label_repel = T, label.size = 6, label_repulsion = 80, label_segment_color = "red", label_point_color = "black",
                    show_stat = T, cells.highlight = T, theme_use = "theme_blank", legend.position = "none", combine = F)
pdf("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/each.timepoint.umap.pdf", width = 8, height = 8)
walk(plots, ~print(.))
dev.off()

plots = CellDimPlot(srat.merge, group.by = "full.name", split.by = "timepoint", reduction = "UMAP", palcolor = colors.lung.merge[levels(srat.merge$split_ctp)],
                    label = F, label_insitu = T, label_repel = T, label.size = 6, label_repulsion = 20, label_segment_color = "red", label_point_color = "black",
                    show_stat = T, cells.highlight = T, theme_use = "theme_blank", legend.position = "top", legend.direction = "horizontal", combine = F)
pdf("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/each.timepoint.full.name.umap.pdf", width = 8*1.5, height = 8*1.5)
walk(plots, ~print(.))
dev.off()

write.csv(srat.merge@meta.data, file = "/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/cache-merge_annotation/metadata.csv")
# run adjust.pictures.R
saveRDS(srat.merge, "/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/cache-merge_annotation/srat.merge.annotated.rds")
# srat.merge = readRDS("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/cache-merge_annotation/srat.merge.annotated.rds")
###############################################################################################################################################
srat.merge = SetIdent(srat.merge, value = "split_ctp")
markers.split_ctp = FindAllMarkers(srat.merge, logfc.threshold = 0.15, test.use = "wilcox", slot = "data", min.pct = 0.15, only.pos = T)
write.csv(markers.split_ctp, file = "/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/cache-merge_annotation/DEGs.csv")
# markers.split_ctp = read.csv("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/cache-merge_annotation/DEGs.csv", row.names = 1)

# 热图修改？
# markers.split_ctp.valid = unique(markers.split_ctp$gene[markers.split_ctp$avg_log2FC > 0.8 & markers.split_ctp$p_val_adj < 0.01])
# markers.split_ctp.up.top = markers.split_ctp %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
# pdf(paste(anno_path, "DEGs.up.pdf", sep = "/"), width = 8*3, height = 8*2)
# Doheatmap_genger(srat.merge, assay = "SCT", slot = "data", sort.by = "split_ctp", split.by = "split_ctp", scale.split.by = NULL, marker.genes = unique(markers.split_ctp.valid),
#                  markers.highlight = unique(markers.split_ctp.up.top$gene), highlight.line.color = "grey", highlight.line.lwd = 0.3, font.size.row = 5, col.name.rot = 90, col.name.size = 8,
#                  annotation.vars = c("timepoint", "nCount_Spatial"), annotation.cols = list("timepoint" = colors.timepoint, "nCount_Spatial" = c("blue", "white", "red")),
#                  heat.min = NULL, heat.max = NULL, heat.scale = seq(-2, 2, length = 100), heat.cols = viridis(n = 100, option = "D"), plot.mean = F)
# dev.off()

###############################################################################################################################################
plots1 = map(meta.info$name, ~SpatialDimPlot_new(srat.merge, group.by = "split_ctp", image = ., image.alpha = 1, crop = F,
                                                 cols = unname(colors.lung.merge[levels(droplevels(srat.merge$split_ctp[srat.merge$orig.ident == .]))]),
                                                 pt.size.factor = 0.1, spot.alpha = 0.8, label = F, show.sb = T, dark = F, draw_channel = F))
plots2 = map(meta.info$name, ~SpatialDimPlot_new(srat.merge, group.by = "split_ctp", image = ., image.alpha = 0, crop = F,
                                                 cols = unname(colors.lung.merge[levels(droplevels(srat.merge$split_ctp[srat.merge$orig.ident == .]))]),
                                                 pt.size.factor = 0.17, label = F, show.sb = F, dark = F, draw_channel = F))
pdf("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/annotation.slice.pdf", width = 8*1.2, height = 8)
walk(plots1, ~print(.))
dev.off()

pdf("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/annotation.cartoon.pdf", width = 8*1.2, height = 8)
walk(plots2, ~print(.))
dev.off()

tmp_path = "/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/light.each.celltype"
for (clus in levels(srat.merge$split_ctp)) {
  message(clus)
  
  srat.merge$tmp = as.character(srat.merge$split_ctp)
  srat.merge$tmp[srat.merge$split_ctp != clus] = "other"
  srat.merge$tmp = factor(srat.merge$tmp, levels = c(clus, "other"))
  plots1 = map(meta.info$name, ~SpatialDimPlot_new(srat.merge, group.by = "tmp", image = ., image.alpha = 0, crop = F,
                                                   cols = unname(colors.lung.merge[levels(srat.merge$tmp)]), spot.alpha = 1,
                                                   pt.size.factor = 0.17, label = F, show.sb = F, dark = F, draw_channel = F))
  srat.tmp = SubsetSTData(srat.merge, idents = clus)
  plots2 = map(meta.info$name, ~SpatialDimPlot_new(srat.tmp, group.by = "split_ctp", image = ., image.alpha = 1, crop = F,
                                                   cols = unname(colors.lung.merge[clus]), spot.alpha = 0.8,
                                                   pt.size.factor = 0.1, label = F, show.sb = F, dark = F, draw_channel = F))
  CheckGC(); gc()

  pdf(paste(tmp_path, paste0(clus, ".cartoon.pdf"), sep = "/"), width = 8*1.2, height = 8)
  walk(plots1, ~print(.))
  dev.off()
  
  pdf(paste(tmp_path, paste0(clus, ".slice.pdf"), sep = "/"), width = 8*1.2, height = 8)
  walk(plots2, ~print(.))
  dev.off()
}

genes.light = unique(unlist(gene.sets.up))
tmp_path = "/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/figure-merge_annotation/light.classic.gene"
for (gene in genes.light) {
  if (gene %in% rownames(srat.merge)) {
    message(gene)
    tmp_path2 = paste(tmp_path, paste0(gene, ".SCT.cartoon.pdf"), sep = "/")
    if (!file.exists(tmp_path2)) {
      p1 = map(levels(srat.merge$orig.ident), ~SpatialFeaturePlot_new(srat.merge, assay = "SCT", feature = gene, image = .,
                                                                      image.alpha = 0, cols = viridis::inferno(100), pt.size.factor = 0.17, spot.alpha = c(1, 1), show.sb = F, dark = F))
      p2 = map(levels(srat.merge$orig.ident), ~SpatialFeaturePlot_new(srat.merge, assay = "SCT", feature = gene, image = .,
                                                                      image.alpha = 1, cols = viridis::inferno(100), pt.size.factor = 0.1, spot.alpha = c(0, 1), show.sb = F, dark = F))
      p3 = map(levels(srat.merge$orig.ident), ~SpatialFeaturePlot_new(srat.merge, assay = "MAGIC_SCT", feature = gene, image = ., 
                                                                      image.alpha = 0, cols = viridis::inferno(100), pt.size.factor = 0.17, spot.alpha = c(1, 1), show.sb = F, dark = F))
      p4 = map(levels(srat.merge$orig.ident), ~SpatialFeaturePlot_new(srat.merge, assay = "MAGIC_SCT", feature = gene, image = .,
                                                                      image.alpha = 1, cols = viridis::inferno(100), pt.size.factor = 0.1, spot.alpha = c(0, 1), show.sb = F, dark = F))
      pdf(tmp_path2, width = 8*1.2, height = 8)
      walk(p1, ~print(.))
      dev.off()
    } else {
      message(paste(gene, "already output!"))
    }
    
    tmp_path2 = paste(tmp_path, paste0(gene, ".SCT.slice.pdf"), sep = "/")
    if (!file.exists(tmp_path2)) {
      pdf(tmp_path2, width = 8*1.2, height = 8)
      walk(p2, ~print(.))
      dev.off()
    }
    
    tmp_path2 = paste(tmp_path, paste0(gene, ".MAGIC.cartoon.pdf"), sep = "/")
    if (!file.exists(tmp_path2)) {
      pdf(tmp_path2, width = 8*1.2, height = 8)
      walk(p3, ~print(.))
      dev.off()
    }
    
    tmp_path2 = paste(tmp_path, paste0(gene, ".MAGIC.slice.pdf"), sep = "/")
    if (!file.exists(tmp_path2)) {
      pdf(tmp_path2, width = 8*1.2, height = 8)
      walk(p4, ~print(.))
      dev.off()
    }
    
  } else {
    message(paste(gene, "does not exists in the assay!"))
  }
}

### output .h5ad ###
seurat2scanpy(srat.merge, out_path = "/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/h5ad", scanpy_name = "srat.merge.annotated.h5ad")

for (i in str_sort(unique(srat.merge$orig.ident), numeric = T)) {
  message(i)
  print(sum(rowSums(srat.merge@assays$Spatial@counts[ , colnames(srat.merge)[srat.merge$orig.ident== i]]) > 0))
}

