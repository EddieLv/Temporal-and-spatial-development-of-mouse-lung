library(getopt)
command=matrix(c("sample","a",1,"character",
                 "slice","b",1,"character",
                 "square_file","c",1,"character",
                 "step","s",1,"character",
                 "mt","m",1,"character",
                 "ribo","r",1,"character",
                 "hb","x",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$sample) || is.null(args$slice) || is.null(args$square_file)|| is.null(args$step)|| is.null(args$mt)|| is.null(args$ribo)|| is.null(args$hb)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

source("/media/biogenger/D/scripts/Spatial-Transcriptome/R/initialize.genger_env.R")
initialize.genger_env()

tmp_func = function(seurat, nCount.high, nCount.low, nFeature.high, nFeature.low, mt.high, ribo.high, hb.high, 
                    nCount_max, nCount_min, nFeature_max, nFeature_min, pct.mt_max, pct.ribo_max, pct.hb_max) {
  params = c(nCount.high, nCount.low, nFeature.high, nFeature.low, mt.high, ribo.high, hb.high)
  names(params) = c("nCount.high", "nCount.low", "nFeature.high", "nFeature.low", "mt.high", "ribo.high", "hb.high")
  label = paste(names(params)[params], collapse = "-")
  
  if (nCount.high) {
    cells1 = seurat$nCount_Spatial >= nCount_max
  } else {
    cells1 = seurat$nCount_Spatial < nCount_max
  }
  
  if (nCount.low) {
    cells2 = seurat$nCount_Spatial < nCount_min
  } else {
    cells2 = seurat$nCount_Spatial >= nCount_min
  }
  
  if (nFeature.high) {
    cells3 = seurat$nFeature_Spatial >= nFeature_max
  } else {
    cells3 = seurat$nFeature_Spatial < nFeature_max
  }
  
  if (nFeature.low) {
    cells4 = seurat$nFeature_Spatial < nFeature_min
  } else {
    cells4 = seurat$nFeature_Spatial >= nFeature_min
  }
  
  if (mt.high) {
    cells5 = seurat$pct.mt >= pct.mt_max
  } else {
    cells5 = seurat$pct.mt < pct.mt_max
  }
  
  if (ribo.high) {
    cells6 = seurat$pct.ribo >= pct.ribo_max
  } else {
    cells6 = seurat$pct.ribo < pct.ribo_max
  }
  
  if (hb.high) {
    cells7 = seurat$pct.hb >= pct.hb_max
  } else {
    cells7 = seurat$pct.hb < pct.hb_max
  }
  
  seurat@meta.data[rownames(seurat@meta.data)[cells1 & cells2 & cells3 & cells4 & cells5 & cells6 & cells7], "labels"] = label
  
  return(seurat)
}

### 0.Set the parameter ###
# Basic path
sample = args$sample
name = args$slice
project_path = "/media/biogenger/D/Projects/CMY"
script_path = "/media/biogenger/D/scripts/Spatial-Transcriptome"
# Upstream path
barcode_path = paste(project_path, "barcodes", sample, "barcodes_A.txt", sep = "/")
countmat_path = paste(project_path, "pipeline_res", sample, "4_zUMIs/zUMIs_output/expression", paste0(sample, ".symbol.tsv"), sep = "/")
spatial_path = paste(project_path, "images", sample, sep = "/")
# Downstream path
qc_path = paste(project_path, "rawQC_res", sample, sep = "/")

### 1.Load data & Create Seurat.object ###
# Countmat file
count.m = t(as.matrix(read.table(countmat_path, check.names = F))) # This is important, or your feature name may be replaced by (-) -> (.)
# Seurat.object
min.cell = ifelse(dim(count.m)[2] / 1000 > 3, round(dim(count.m)[2] / 1000), 3)
srat = CreateSeuratObject(counts = count.m, project = name, assay = "Spatial", min.cells = min.cell, min.features = 0)
print(paste("min.cell is", min.cell))

### 2.Add spatial.info ###
# Create spatial/tissue_positions_list.csv
tissue_positions_list = paste(spatial_path, "tissue_positions_list.csv", sep = "/")
if (! file.exists(tissue_positions_list)) {
  print(paste(tissue_positions_list, "does not exist, it will be created!"))
  create_tissue_position(ifelse(grepl("E12.5", sample) == 1, 50, 96), tissue_positions_list)
}
# Add spatial.info
image = Read10X_Image(spatial_path, image.name = paste0("FIX&", sample, ".png"))
image = image[Cells(srat)]
DefaultAssay(image) = "Spatial"
srat[[name]] = image
srat = AddMetaData(srat, metadata = subset(srat@images[[name]]@coordinates, select = c("row", "col", "imagerow", "imagecol")), 
                   col.name = c("barcodeB", "barcodeA", "imagerow", "imagecol"))
srat[["sample"]] = sample
srat[["timepoint"]] = sapply(str_split(srat$orig.ident, "_"), function(x){x[1]})
srat[["slice"]] = sapply(str_split(srat$orig.ident, "_"), function(x){x[2]})
# Add STutility.info
###########
info.table = data.frame(samples = c(countmat_path), 
                        imgs = c(paste(spatial_path, paste0("FIX&", sample, ".png"), sep = "/")), 
                        spotfiles = c(paste(spatial_path, "tissue_positions_list.csv", sep = "/")), 
                        json = c(paste(spatial_path, "scalefactors_json.json", sep = "/")), 
                        tech = c("DBIT-seq"))
se = InputFromTable(info.table, 
                    minUMICountsPerGene = 0, 
                    minSpotsPerGene = min.cell, 
                    minGenesPerSpot = 0, 
                    minUMICountsPerSpot = 0, 
                    platform = "Visium")
srat@tools$Staffli = se@tools$Staffli
rownames(srat@tools$Staffli@meta.data) = sapply(str_split(rownames(srat@tools$Staffli@meta.data), "_"), function(x){x[1]})
srat@tools$Staffli@meta.data$id = 1:dim(srat@tools$Staffli@meta.data)[1]
srat@tools$Staffli@meta.data = srat@tools$Staffli@meta.data[rownames(srat@meta.data), ]
srat[["id"]] = srat@tools$Staffli@meta.data[rownames(srat@meta.data), "id"]
srat[["labels"]] = "Default"
srat = LoadImages(srat, time.resolve = F, verbose = TRUE)

### 3.Quality Control ###
if (! dir.exists(qc_path)) dir.create(qc_path, recursive = T)
# Remove blank spots
manual.remove.spot = read.table(paste(spatial_path, "manual.remove.txt", sep = "/"))
srat.sub = SubsetSTData(srat, spots = manual.remove.spot$V1, invert = T)
square.remove.spot = read.table(args$square_file)
srat.sub = SubsetSTData(srat.sub, spots = square.remove.spot$V1, invert = T)
# Other variable
srat.sub[["pct.mt"]] = PercentageFeatureSet(srat.sub, pattern = "^mt-")
srat.sub[["pct.ribo"]] = PercentageFeatureSet(srat.sub, pattern = "^Rp[sl]")
srat.sub[["pct.hb"]] = PercentageFeatureSet(srat.sub, pattern = "^Hb[ab]-")
# Log-normalize nCount & nFeature
srat.sub@meta.data = srat.sub@meta.data %>% mutate(nCount_log10 = log10(nCount_Spatial), nFeature_log10 = log10(nFeature_Spatial))
p1.1 = SpatialFeaturePlot_new(srat.sub, feature = "nCount_Spatial", image = name, image.alpha = 0, spot.alpha = c(1, 1), 
                              pt.size.factor = 0.17, draw_channel = F, show.sb = T, dark = T)
p1.2 = SpatialFeaturePlot_new(srat.sub, feature = "nCount_log10", image = name, image.alpha = 0, spot.alpha = c(1, 1), 
                              pt.size.factor = 0.17, draw_channel = F, show.sb = T, dark = T)
p2.1 = SpatialFeaturePlot_new(srat.sub, feature = "nCount_Spatial", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                              draw_channel = F, show.sb = T, dark = T)
p2.2 = SpatialFeaturePlot_new(srat.sub, feature = "nCount_log10", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                              draw_channel = F, show.sb = T, dark = T)
p3.1 = SpatialFeaturePlot_new(srat.sub, feature = "nFeature_Spatial", image = name, image.alpha = 0, spot.alpha = c(1, 1), 
                              pt.size.factor = 0.17, draw_channel = F, show.sb = T, dark = T)
p3.2 = SpatialFeaturePlot_new(srat.sub, feature = "nFeature_log10", image = name, image.alpha = 0, spot.alpha = c(1, 1), 
                              pt.size.factor = 0.17, draw_channel = F, show.sb = T, dark = T)
p4.1 = SpatialFeaturePlot_new(srat.sub, feature = "nFeature_Spatial", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                              draw_channel = F, show.sb = T, dark = T)
p4.2 = SpatialFeaturePlot_new(srat.sub, feature = "nFeature_log10", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                              draw_channel = F, show.sb = T, dark = T)
p5.1 = SpatialFeaturePlot_new(srat.sub, feature = "pct.mt", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                              draw_channel = F, show.sb = T, dark = T)
p5.2 = SpatialFeaturePlot_new(srat.sub, feature = "pct.ribo", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                              draw_channel = F, show.sb = T, dark = T)
p5.3 = SpatialFeaturePlot_new(srat.sub, feature = "pct.hb", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                              draw_channel = F, show.sb = T, dark = T)
pdf(paste(qc_path, "rawQC.pdf", sep = "/"), width = 18, height = 12)
plot_QC(srat.sub)
dev.off()
pdf(paste(qc_path, "rawQC_slice.pdf", sep = "/"), width = 8*2, height = 8)
wrap_plots(p1.1, p1.2)
wrap_plots(p2.1, p2.2)
wrap_plots(p3.1, p3.2)
wrap_plots(p4.1, p4.2)
wrap_plots(p5.1, p5.2, p5.3)
dev.off()

# Filter Again
# nCount_min = quantile(srat.sub[["nCount_Spatial"]][[1]], probs = 0.005)[[1]]
# nCount_max = quantile(srat.sub[["nCount_Spatial"]][[1]], probs = 0.995)[[1]]
# nFeature_min = quantile(srat.sub[["nFeature_Spatial"]][[1]], probs = 0.005)[[1]]
# nFeature_max = quantile(srat.sub[["nFeature_Spatial"]][[1]], probs = 0.995)[[1]]
# 每个细胞平均检测到500个gene具有生物学意义，不依赖于数据的分布！
nCount_min = 1000
nCount_max = quantile(srat.sub[["nCount_Spatial"]][[1]], probs = 0.99)[[1]]
nFeature_min = 500
nFeature_max = quantile(srat.sub[["nFeature_Spatial"]][[1]], probs = 0.99)[[1]]
pct.mt_max = as.numeric(args$mt); pct.ribo_max = as.numeric(args$ribo); pct.hb_max = as.numeric(args$hb)
print(paste("nCount_max", nCount_max))
print(paste("nFeature_max", nFeature_max))
print(paste("pct.mt_max", pct.mt_max))
print(paste("pct.ribo_max", pct.ribo_max))
print(paste("pct.hb_max", pct.hb_max))
# Only one
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
# Both
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = T, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
# Tripple
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = T, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = T, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = T, mt.high = T, ribo.high = F, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = T, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = T, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = T, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
# Four
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = T, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = T, nFeature.low = F, mt.high = T, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = T, nFeature.low = F, mt.high = T, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = T, nFeature.low = F, mt.high = F, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = T, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = T, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = T, nCount.low = F, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = T, mt.high = T, ribo.high = T, hb.high = F, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = T, mt.high = T, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)
srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = T, nFeature.high = F, nFeature.low = T, mt.high = F, ribo.high = F, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

srat.sub = tmp_func(srat.sub, nCount.high = F, nCount.low = F, nFeature.high = F, nFeature.low = F, mt.high = T, ribo.high = T, hb.high = T, 
                    nCount_max = nCount_max, nCount_min = nCount_min, nFeature_max = nFeature_max, nFeature_min = nFeature_min, pct.mt_max = pct.mt_max, pct.ribo_max = pct.ribo_max, pct.hb_max = pct.hb_max)

table(srat.sub$labels)
srat.sub.cutoff = srat.sub[ , nCount_min < srat.sub$nCount_Spatial & srat.sub$nCount_Spatial < nCount_max & 
                            nFeature_min < srat.sub$nFeature_Spatial & srat.sub$nFeature_Spatial < nFeature_max & 
                            srat.sub$pct.mt < pct.mt_max & 
                            srat.sub$pct.ribo < pct.ribo_max & 
                            srat.sub$pct.hb < pct.hb_max]
# Crop manually if needed
#srat.sub = shiny_crop_gener(srat.sub, type = NULL, res = 1080, verbose = F)
p = SpatialDimPlot_new(srat.sub, group.by = "labels", image = name, image.alpha = 1, crop = F, show.sb = T, dark = F)
pdf(paste(qc_path, "mt&hb.situation.pdf", sep = "/"), width = 8, height = 8)
print(p)
dev.off()

srat.sub.cutoff[["pct.mt"]] = PercentageFeatureSet(srat.sub.cutoff, pattern = "^mt-")
srat.sub.cutoff[["pct.ribo"]] = PercentageFeatureSet(srat.sub.cutoff, pattern = "^Rp[sl]")
srat.sub.cutoff[["pct.hb"]] = PercentageFeatureSet(srat.sub.cutoff, pattern = "^Hb[ab]-")
pdf(paste(qc_path, "cutoffQC.pdf", sep = "/"), width = 18, height = 12)
plot_QC(srat.sub.cutoff)
dev.off()
  
if (args$step == "whole") {
  ### 4.Find doublet ###
  methods = c("cxds", "bcds", "hybrid", "scDblFinder", "Scrublet", "DoubletDetection", "DoubletFinder")
  srat.sub.doublet = find_doublet(srat.sub.cutoff, methods = methods, outPath = paste(qc_path, "find_Doublet_multi.svg", sep = "/"), step = "all")
  # 全排列三三交叉Doublet
  pd_doublets = rownames(srat.sub.doublet@meta.data)[unique(c(Reduce(intersect, list(which(srat.sub.doublet@meta.data$pd_scDblFinder == "doublet"), 
                                                                                     which(srat.sub.doublet@meta.data$pd_DoubletDetection == "doublet"), 
                                                                                     which(srat.sub.doublet@meta.data$pd_DoubletFinder == "doublet"))), 
                                                              Reduce(intersect, list(which(srat.sub.doublet@meta.data$pd_scDblFinder == "doublet"), 
                                                                                     which(srat.sub.doublet@meta.data$pd_DoubletDetection == "doublet"), 
                                                                                     which(srat.sub.doublet@meta.data$Scrublet == "doublet"))), 
                                                              Reduce(intersect, list(which(srat.sub.doublet@meta.data$pd_scDblFinder == "doublet"), 
                                                                                     which(srat.sub.doublet@meta.data$pd_DoubletFinder == "doublet"), 
                                                                                     which(srat.sub.doublet@meta.data$Scrublet == "doublet"))), 
                                                              Reduce(intersect, list(which(srat.sub.doublet@meta.data$pd_DoubletDetection == "doublet"), 
                                                                                     which(srat.sub.doublet@meta.data$pd_DoubletFinder == "doublet"), 
                                                                                     which(srat.sub.doublet@meta.data$Scrublet == "doublet")))
  ))]
  srat.sub.doublet@meta.data$pd_doublets = "singlet"
  srat.sub.doublet@meta.data$pd_doublets[match(pd_doublets, Cells(srat.sub.doublet))] = "doublet"
  p = SpatialDimPlot_new(srat.sub.doublet, group.by = "pd_doublets", image = name, image.alpha = 1, crop = F, show.sb = T, dark = F)
  pdf(paste(qc_path, "doublet.situation.pdf", sep = "/"), width = 8*2, height = 8)
  print(p)
  dev.off()
  # Filter Doublet
  srat.sub.final = SubsetSTData(srat.sub.doublet, spots = match(pd_doublets, Cells(srat.sub.doublet)), invert = T)
  print("Stat of final")
  print(srat.sub.final)
  print(srat.sub.final@meta.data %>% dplyr::select(nCount_Spatial, nFeature_Spatial, pct.mt, pct.ribo, pct.hb) %>% summary())
  
  ### 5.Remove useless gene ###
  # genes_remove = grep("^Hba-|^Hbb-|^mt-|^Rps|^Rpl", x = rownames(srat.sub.final), value = T)
  genes_remove = gene.table.remove$`gene name`
  print(paste("Removing", length(intersect(rownames(srat.sub.final), genes_remove)), "genes useless!"))
  genes_keep = setdiff(rownames(srat.sub.final), genes_remove)
  srat.sub.final = srat.sub.final[genes_keep, ]
  
  genes_remove_zero_exp = rownames(srat.sub.final)[rowSums(srat.sub.final) == 0]
  print(paste("Removing", length(intersect(rownames(srat.sub.final), genes_remove_zero_exp)), "genes of zero expression"))
  genes_keep = setdiff(rownames(srat.sub.final), genes_remove_zero_exp)
  srat.sub.final = srat.sub.final[genes_keep, ]
  
  spot_remove_zero_exp = colnames(srat.sub.final)[colSums(srat.sub.final) == 0]
  print(paste("Removing", length(intersect(colnames(srat.sub.final), spot_remove_zero_exp)), "spots of zero expression"))
  spots_keep = setdiff(colnames(srat.sub.final), spot_remove_zero_exp)
  srat.sub.final = srat.sub.final[ , spots_keep]
  
  srat.sub.final[["pct.mt"]] = PercentageFeatureSet(srat.sub.final, pattern = "^mt-")
  srat.sub.final[["pct.ribo"]] = PercentageFeatureSet(srat.sub.final, pattern = "^Rp[sl]")
  srat.sub.final[["pct.hb"]] = PercentageFeatureSet(srat.sub.final, pattern = "^Hb[ab]-")
  
  # Log-normalize nCount & nFeature
  srat.sub.final@meta.data = srat.sub.final@meta.data %>% mutate(nCount_log10 = log10(nCount_Spatial), nFeature_log10 = log10(nFeature_Spatial))
  p1.1 = SpatialFeaturePlot_new(srat.sub.final, feature = "nCount_Spatial", image = name, image.alpha = 0, spot.alpha = c(1, 1), 
                                pt.size.factor = 0.17, draw_channel = F, show.sb = T, dark = T)
  p1.2 = SpatialFeaturePlot_new(srat.sub.final, feature = "nCount_log10", image = name, image.alpha = 0, spot.alpha = c(1, 1), 
                                pt.size.factor = 0.17, draw_channel = F, show.sb = T, dark = T)
  p2.1 = SpatialFeaturePlot_new(srat.sub.final, feature = "nCount_Spatial", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                                draw_channel = F, show.sb = T, dark = T)
  p2.2 = SpatialFeaturePlot_new(srat.sub.final, feature = "nCount_log10", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                                draw_channel = F, show.sb = T, dark = T)
  p3.1 = SpatialFeaturePlot_new(srat.sub.final, feature = "nFeature_Spatial", image = name, image.alpha = 0, spot.alpha = c(1, 1), 
                                pt.size.factor = 0.17, draw_channel = F, show.sb = T, dark = T)
  p3.2 = SpatialFeaturePlot_new(srat.sub.final, feature = "nFeature_log10", image = name, image.alpha = 0, spot.alpha = c(1, 1), 
                                pt.size.factor = 0.17, draw_channel = F, show.sb = T, dark = T)
  p4.1 = SpatialFeaturePlot_new(srat.sub.final, feature = "nFeature_Spatial", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                                draw_channel = F, show.sb = T, dark = T)
  p4.2 = SpatialFeaturePlot_new(srat.sub.final, feature = "nFeature_log10", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                                draw_channel = F, show.sb = T, dark = T)
  p5.1 = SpatialFeaturePlot_new(srat.sub.final, feature = "pct.mt", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                                draw_channel = F, show.sb = T, dark = T)
  p5.2 = SpatialFeaturePlot_new(srat.sub.final, feature = "pct.ribo", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                                draw_channel = F, show.sb = T, dark = T)
  p5.3 = SpatialFeaturePlot_new(srat.sub.final, feature = "pct.hb", image = name, image.alpha = 0.8, spot.alpha = c(0, 1), 
                                draw_channel = F, show.sb = T, dark = T)
  pdf(paste(qc_path, "finalQC.pdf", sep = "/"), width = 18, height = 12)
  print(plot_QC(srat.sub.final))
  dev.off()
  pdf(paste(qc_path, "finalQC_slice.pdf", sep = "/"), width = 8*2, height = 8)
  print(wrap_plots(p1.1, p1.2))
  print(wrap_plots(p2.1, p2.2))
  print(wrap_plots(p3.1, p3.2))
  print(wrap_plots(p4.1, p4.2))
  print(wrap_plots(p5.1, p5.2, p5.3))
  dev.off()
  
  print("Stat of remove zero-expressed genes/spots")
  print(srat.sub.final)
  print(srat.sub.final@meta.data %>% dplyr::select(nCount_Spatial, nFeature_Spatial, pct.mt, pct.ribo, pct.hb) %>% summary())
  write.table(t(as.matrix(srat.sub.final@assays$Spatial@counts)), file = paste(qc_path, "final.count.tsv", sep = "/"), 
              row.names = T, sep = "\t", quote = F)
  write.csv(srat.sub.final@meta.data, file = paste(qc_path, "metadata.csv", sep = "/"), row.names = T)
  saveRDS(srat.sub.final, file = paste(qc_path, "srat.raw.qc.rds", sep = "/"))
}

### 6.SoupX ###
#srat.soup = SCTransform(srat.sub.final, assay = "Spatial", vst.flavor = "v2", verbose = F) %>% 
#  RunPCA(assay = "SCT", verbose = F)
#pc = Stdev(srat.soup, reduction = "pca")
#eigValues = pc ** 2
#varExplained = eigValues / sum(eigValues)
#for (n in 1:50) {
#  if (sum(varExplained[1:n]) >= 0.92) {
#    print(n)
#    break
#  }
#}
#srat.soup = RunUMAP(srat.soup, reduction = "pca", dims = 1:n, n.neighbors = 50) %>% 
#  FindNeighbors(reduction = "pca", dims = 1:n) %>% 
#  FindClusters(algorithm = 1, resolution = 0.8)
#meta.data = srat.soup@meta.data
#toc = as.matrix(srat.sub.final@assays$Spatial@counts)
#tod = count.m[rownames(toc), ]
#srat.soup = SoupChannel(tod, toc, calcSoupProfile = F)
#droplet.freq = colSums(count.m)[manual.remove.spot$V1]
#droplet.mean = mean(droplet.freq, na.rm = T)
#droplet.median = median(droplet.freq, na.rm = T)
#droplet.cutoff = min(min(srat.sub.final$nCount_Spatial), 0.5 * droplet.mean + 0.5 * droplet.median)
#print(paste("droplet.mean:", droplet.mean, "droplet.median:", droplet.median, "droplet.cutoff:", droplet.cutoff))
#srat.soup = estimateSoup(srat.soup, soupRange = c(0, droplet.cutoff), keepDroplets = F)
#srat.soup = setClusters(srat.soup, setNames(meta.data$seurat_clusters, rownames(meta.data)))
#rho = max(0.2, autoEstCont(srat.soup)$fit$rhoEst)
#print(paste("selected rho:", rho))
#srat.soup = setContaminationFraction(srat.soup, rho)
## 校正矩阵
#mat.corrected = as.matrix(adjustCounts(srat.soup))
#write.table(t(mat.corrected), file = paste(qc_path, "final.soupX.count.tsv", sep = "/"), 
#            row.names=T, sep = "\t", quote = F)
#save(srat.sub.final, file = paste(qc_path, "srat.soup.Rdata", sep = "/"))




##############################################################################################################################
# 对比SoupX的效果: Epcam(Epith) Ddc(PNEC) Pecam1(Endo) Pdgfra(Immune) Acta2(Mesen)
#nosoupx = SCTransform(srat.sub.final, assay = "Spatial", vst.flavor = "v2", verbose = F) %>%
#  RunPCA(assay = "SCT", verbose = F)
#pc = Stdev(nosoupx, reduction = "pca")
#eigValues = pc ** 2
#varExplained = eigValues / sum(eigValues)
#for (n in 1:50) {
#  if (sum(varExplained[1:n]) >= 0.92) {
#    print(n)
#    break
#  }
#}
#nosoupx = RunUMAP(nosoupx, reduction = "pca", dims = 1:n, n.neighbors = 50) %>%
#  FindNeighbors(reduction = "pca", dims = 1:n) %>%
#  FindClusters(algorithm = 1, resolution = 0.8)


#soupx = CreateSeuratObject(counts = mat.corrected, project = name, assay = "Spatial", min.cells = 0, min.features = 0) %>%
#  SCTransform(assay = "Spatial", vst.flavor = "v2", verbose = F) %>%
#  RunPCA(assay = "SCT", verbose = F)
#pc = Stdev(soupx, reduction = "pca")
#eigValues = pc ** 2
#varExplained = eigValues / sum(eigValues)
#for (n in 1:50) {
#  if (sum(varExplained[1:n]) >= 0.92) {
#    print(n)
#    break
#  }
#}
#soupx = RunUMAP(soupx, reduction = "pca", dims = 1:n, n.neighbors = 50) %>%
#  FindNeighbors(reduction = "pca", dims = 1:n) %>%
#  FindClusters(algorithm = 1, resolution = 0.8)

#p1 = VlnPlot(nosoupx, features = c("Epcam", "Pecam1", "Acta2", "Ptprc"), ncol = 4, slot = "counts")
#p2 = VlnPlot(soupx, features = c("Epcam", "Pecam1", "Acta2", "Ptprc"), ncol = 4, slot = "counts")
#pdf(paste(qc_path, "SoupX.pdf", sep = "/"), width = 8*4, height = 8)
#p1
#p2
#dev.off()

