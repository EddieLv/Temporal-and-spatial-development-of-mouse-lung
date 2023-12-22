library(getopt)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
command=matrix(c("pixel_stat_file","f",1,"character",
                 "work_path","p",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$pixel_stat_file) || is.null(args$work_path)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

work_path = args$work_path # "/media/biogenger/D/Projects/CZP/Cere-24-20220416/"
pixel_stat_file = args$pixel_stat_file # "/media/biogenger/D/Projects/CZP/Cere-24-20220416/map_pixel.csv"

setwd(work_path)

stat = read.csv(pixel_stat_file)
stat=stat[,c("iB","iA","raw_count")]
nchannel = max(stat$iB)

if (nchannel == 50) {
  all_combi.m=spread(stat,iB,raw_count) # Attention: makes the matrix row belongs to bA(1->50); makes the matrix column belongs to bB(1->50)
  all_combi.m=all_combi.m[,-1]
  all_combi.m=as.matrix(all_combi.m[,rev(1:nchannel)]) # reverse the columns of matrix
}else if (nchannel == 96) {
  all_combi.m=spread(stat,iB,raw_count) # Attention: makes the matrix row belongs to bA(1->96); makes the matrix column belongs to bB(1->96)
  all_combi.m=all_combi.m[,-1]
  all_combi.m=as.matrix(all_combi.m[rev(1:nchannel),]) # reverse the rows of matrix
}else {
  print("Error in barcode num!")
}
library(circlize)
svg("original.svg",width = 8,height = 8)
Heatmap(all_combi.m,cluster_rows=F,cluster_columns=F,rect_gp = gpar(col = "white", lwd = 5),
        col = colorRamp2(c(0, max(all_combi.m)/5, max(all_combi.m)), c("blue", "green", "red")),
        show_heatmap_legend = F,show_row_names=F,show_column_names=F)
dev.off()



