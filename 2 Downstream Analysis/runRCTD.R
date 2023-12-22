set.seed(666)
source("/media/biogenger/D/scripts/Spatial-Transcriptome/R/initialize.genger_env.R")
initialize.genger_env()

tmp1 = function(x, num) {
  p = map(num, 
          ~stack.barplot(srat.sub@meta.data, 
                         x.var = paste0("rctd.", x, ".", .*100, ".type"), 
                         y.var = paste0("rctd.", x, ".sig.count.type"), 
                         fill.col = ggsci::pal_aaas()(4), 
                         type = "count") + ggtitle(paste0(x, "_cutoff", .)))
  
  return(wrap_plots(p, ncol = 3))
}

tmp2 = function(seurat, mode, num, slice) {
  if (sum(is.na(seurat[[paste("rctd", mode, num, "weight", sep = ".")]])) == dim(seurat)[2]) {
    message("No valid cells!")
    return(ggplot())
  }
  
  p1 = SpatialDimPlot_new(seurat, group.by = paste("rctd", mode, num, "sig", sep = "."), image = slice)
  p2 = SpatialDimPlot_new(subset(seurat, cells = rownames(seurat@meta.data)[seurat@meta.data[, paste("rctd", mode, num, "sig", sep = ".")] %in% c("conf", "not")]), 
                          group.by = paste("rctd", mode, num, "type", sep = "."), image = slice) + theme(legend.spacing.y = unit(0, "mm")) + guides(fill = guide_legend(ncol = 1, byrow = T))
  if (sum(seurat@meta.data[, paste("rctd", mode, num, "sig", sep = ".")] %in% c("conf")) == 0) {
    message("No valid cells!")
    p3 = ggplot()
  } else {
    p3 = SpatialDimPlot_new(subset(seurat, cells = rownames(seurat@meta.data)[seurat@meta.data[, paste("rctd", mode, num, "sig", sep = ".")] %in% c("conf")]), 
                            group.by = paste("rctd", mode, num, "type", sep = "."), image = slice) + theme(legend.spacing.y = unit(0, "mm")) + guides(fill = guide_legend(ncol = 1, byrow = T))
  }
  
  return(wrap_plots(p1, p2, p3, ncol = 3))
}

srat.merge = readRDS("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/3-merge_annotation/cache-merge_annotation/srat.merge.annotated.rds")
# srat.merge = readRDS("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig1/overall_umap/10um/2-merge_annotation/cache-merge_annotation/srat.merge.annotated.rds")

### 0.Initialize ###
sample.list = c("E12.5_slice13", "E12.5_slice21", "E15.5_slice7", "E16.5_slice24", "E16.5_slice39", "E17.5_slice13", "E18.5_slice4", "E18.5_slice12", "P0_slice12", "P3_slice1")
# sample.list = c("E16.5_slice18", "E18.5_slice31")
reference.list = c("E12", "E12", "E15", "E16", "E16", "E18", "E18", "E18", "P0", "P3")
# reference.list = c("E16", "E18")
size.list = c(0.5, 0.5, 0.3, 0.7, 0.7, 0.5, 0.5, 0.5, 0.5)
# size.list = c(0.4, 0.4)

library(spacexr)
for (i in 1:length(sample.list)) {
  sample = sample.list[i]
  print(sample)
  reference = reference.list[i]
  ################################################################################################################################
  ### Deconvolution ###
  rctd_path = "/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/Overall/Fig2/Deconvolution/RCTD/result"
  advance_path = paste(rctd_path, sample, sep = "/")
  if (! dir.exists(advance_path)) dir.create(advance_path, recursive = T)
  # Load SingleCell Reference
  srat.sc = readRDS(paste("/media/biogenger/D/Projects/CMY/Analysis/mouse_lung/scRNA_reference/cache-split_time_reannotation", paste0(reference, ".rds"), sep = "/"))
  srat.sc$split_reanno_ctp = factor(as.character(srat.sc$split_reanno_ctp), levels = unique(srat.sc$split_reanno_ctp))
  sc.ref = Reference(srat.sc@assays$RNA@counts, cell_types = setNames(as.factor(srat.sc$split_reanno_ctp), nm = colnames(srat.sc)), require_int = F, min_UMI = 0)

  saveRDS(sc.ref, file = paste(advance_path, "sc.ref.rctd.rds", sep = "/"))
  #sc.ref = readRDS(paste(advance_path, "sc.ref.rctd.rds", sep = "/"))
  
  # Load ST Data
  srat.sub = SubsetSTData(srat.merge, spots = rownames(srat.merge@meta.data)[srat.merge$orig.ident == sample])
  puck = SpatialRNA(counts = srat.sub@assays$Spatial@counts, 
                    coords = GetTissueCoordinates(srat.sub, image = sample) %>% set_names("x", "y"), 
                    nUMI = srat.sub$nCount_Spatial)
  
  myRCTD.reps = create.RCTD(puck, sc.ref, CELL_MIN_INSTANCE = 1, max_cores = 12)
  myRCTD.reps.doublet = run.RCTD(myRCTD.reps, doublet_mode = "doublet")
  myRCTD.reps.full = run.RCTD(myRCTD.reps, doublet_mode = "full")
  
  # @ 3 modes
  # ‘full mode’ (no restrictions on number of cell types)
  # ‘multi mode’ (finitely many cell types per pixel, e.g. 3 or 4).
  # ‘doublet mode’ (at most 1-2 cell types per pixel)
  
  # @full
  # --@results
  # ----@weights
  
  # @multi
  # --@results
  # ----@[[1]] (pixel)
  # ----@[[2]] (pixel)
  # ----@[[...]] (pixel)
  # ----@[[n]] (pixel)
  
  # @doublet
  # --@results
  # ----@results_df (“singlet” (1 cell type on pixel), “doublet_certain” (2 cell types on pixel), “doublet_uncertain” (2 cell types on pixel, but only confident of 1), “reject” (no prediction given for pixel).)
  # ----@weights (!!!the same as full weights!!!)
  # ----@weights_doublet
  # ----@score_mat
  
  srat.sub = add.rctd2meta(srat.sub, rctd.full = myRCTD.reps.full, rctd.doublet = myRCTD.reps.doublet)
  srat.sub = AddMetaData(srat.sub, metadata = add.count.type(srat.sub@meta.data, levels(sc.ref@cell_types), mode = "full"), col.name = "rctd.full.sig.count.type")
  srat.sub = AddMetaData(srat.sub, metadata = add.count.type(srat.sub@meta.data, levels(sc.ref@cell_types), mode = "doublet"), col.name = "rctd.doublet.sig.count.type")
  
  for (w in seq(0.1, 0.9 ,0.1)) {
    srat.sub = AddMetaData(srat.sub, metadata = add.best.celltype(srat.sub@meta.data, levels(sc.ref@cell_types), mode = "full", max.n = 1, min.weight = w, type = "celltype"), col.name = paste0("rctd.full.", w*100, ".type"))
    srat.sub = AddMetaData(srat.sub, metadata = add.best.celltype(srat.sub@meta.data, levels(sc.ref@cell_types), mode = "full", max.n = 1, min.weight = w, type = "sig"), col.name = paste0("rctd.full.", w*100, ".sig"))
    srat.sub = AddMetaData(srat.sub, metadata = add.best.celltype(srat.sub@meta.data, levels(sc.ref@cell_types), mode = "full", max.n = 1, min.weight = w, type = "weight"), col.name = paste0("rctd.full.", w*100, ".weight"))
    
    srat.sub = AddMetaData(srat.sub, metadata = add.best.celltype(srat.sub@meta.data, levels(sc.ref@cell_types), mode = "doublet", max.n = 1, min.weight = w, type = "celltype"), col.name = paste0("rctd.doublet.", w*100, ".type"))
    srat.sub = AddMetaData(srat.sub, metadata = add.best.celltype(srat.sub@meta.data, levels(sc.ref@cell_types), mode = "doublet", max.n = 1, min.weight = w, type = "sig"), col.name = paste0("rctd.doublet.", w*100, ".sig"))
    srat.sub = AddMetaData(srat.sub, metadata = add.best.celltype(srat.sub@meta.data, levels(sc.ref@cell_types), mode = "doublet", max.n = 1, min.weight = w, type = "weight"), col.name = paste0("rctd.doublet.", w*100, ".weight"))
  }

  
  # @feature
  # paste0("rctd.", feature, ".weight.full")
  # paste0("rctd.", feature, ".sig.multi")
  # paste0("rctd.", feature, ".weight.multi")
  # paste0("rctd.", feature, ".sig.doublet")
  # paste0("rctd.", feature, ".weight.doublet")
  
  saveRDS(myRCTD.reps.full, file = paste(advance_path, "rctd.full.rds", sep = "/"))
  saveRDS(myRCTD.reps.doublet, file = paste(advance_path, "rctd.doublet.rds", sep = "/"))
  saveRDS(srat.sub, file = paste(advance_path, "srat.rctd.rds", sep = "/"))
  
  #myRCTD.reps.full = readRDS(paste(advance_path, "rctd.full.rds", sep = "/"))
  #myRCTD.reps.doublet = readRDS(paste(advance_path, "rctd.doublet.rds", sep = "/"))
  #srat.sub = readRDS(paste(advance_path, "srat.rctd.rds", sep = "/"))
  
  fig_path = paste(advance_path, "plot", sep = "/")
  if (! dir.exists(fig_path)) dir.create(fig_path)
  cell_types = gsub("\\+", ".", levels(sc.ref@cell_types))
  cell_types = gsub(" ", ".", cell_types)
  
  pdf(paste(fig_path, paste0("rctd.celltype.full.", sample, ".pdf"), sep = "/"), width = 8*2, height = 8)
  walk(cell_types, ~print(SpatialPlot.rctd(srat.sub, image = sample, celltype = ., mode = "full", crop = T)))
  dev.off()
  pdf(paste(fig_path, paste0("rctd.celltype.doublet.", sample, ".pdf"), sep = "/"), width = 8*2, height = 8)
  walk(cell_types, ~print(SpatialPlot.rctd(srat.sub, image = sample, celltype = ., mode = "doublet", crop = T)))
  dev.off()
  
  pdf(paste(fig_path, paste0("rctd.pieplot.full.", sample, ".pdf"), sep = "/"), width = 8*1.5, height = 8)
  print(rctd.pieplot(srat.sub, features = paste0("rctd.", cell_types, ".weight.full"), pie.size = size.list[i], show.image = T) + scale_fill_manual(values = pal))
  print(rctd.pieplot(srat.sub, features = paste0("rctd.", cell_types, ".weight.full"), pie.size = size.list[i], show.image = F) + scale_fill_manual(values = pal))
  dev.off()
  pdf(paste(fig_path, paste0("rctd.pieplot.doublet.", sample, ".pdf"), sep = "/"), width = 8*1.5, height = 8)
  print(rctd.pieplot(srat.sub, features = paste0("rctd.", cell_types, ".weight.doublet"), pie.size = size.list[i], show.image = T) + scale_fill_manual(values = pal))
  print(rctd.pieplot(srat.sub, features = paste0("rctd.", cell_types, ".weight.doublet"), pie.size = size.list[i], show.image = F) + scale_fill_manual(values = pal))
  dev.off()
  
  p1 = map(c("full", "doublet"), ~tmp1(., seq(0.1, 0.9, 0.1)))
  pdf(paste(fig_path, paste0("rctd.best.celltype.stack.pdf"), sep = "/"), width = 8*3, height = 8*3)
  print(p1)
  dev.off()
  
  p1 = map(paste0("rctd.full.", seq(10, 90, 10), ".weight"), ~plot_violin(srat.sub, feature = .))
  p2 = map(paste0("rctd.doublet.", seq(10, 90, 10), ".weight"), ~plot_violin(srat.sub, feature = .))
  pdf(paste(fig_path, paste0("rctd.best.celltype.violin.pdf"), sep = "/"), width = 8*3, height = 8*3)
  print(wrap_plots(p1, ncol = 3))
  print(wrap_plots(p2, ncol = 3))
  dev.off()
  
  pdf(paste(fig_path, paste0("rctd.best.celltype.slice.pdf"), sep = "/"), width = 8*3, height = 8)
  walk(seq(10, 90, 10), ~print(tmp2(srat.sub, "full", ., sample)))
  walk(seq(10, 90, 10), ~print(tmp2(srat.sub, "doublet", ., sample)))
  dev.off()
  
  write.csv(srat.sub@meta.data, file = paste(advance_path, paste0("metadata", "_", sample, ".csv"), sep = "/"), row.names = T)
  
  gc()
}
是
