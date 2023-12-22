library(getopt)
library(ggplot2)
library(dplyr)
command=matrix(c("rds_mat","i",1,"character",
                 "barcode_csv","b",1,"character",   
                 "name_mat","o",1,"character",
                 "help","h",0,"logical"),byrow=T,ncol=4)

args=getopt(command)
if (!is.null(args$help) || is.null(args$rds_mat) || is.null(args$barcode_csv) || is.null(args$name_mat)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

rds.mat = args$rds_mat # "/media/biogenger/D/Projects/CMY/pipeline_res/P14_D_3_RNA_cmy_23.1.9/4_zUMIs/zUMIs_output/expression/P14_D_3_RNA_cmy_23.1.9.dgecounts.rds"
barcode.csv = args$barcode_csv # "/media/biogenger/D/Projects/CMY/pipeline_res/P14_D_3_RNA_cmy_23.1.9/3_umi_tools/debarcoded_passed_reads_stat.csv"
name.mat = args$name_mat # "/media/biogenger/D/Projects/CMY/pipeline_res/P14_D_3_RNA_cmy_23.1.9/4_zUMIs/zUMIs_output/expression/P14_D_3_RNA_cmy_23.1.9.symbol.tsv"

gene.table = read.table("/media/biogenger/D/Reference/mouse/GRCm38_mm10/ensembl/Mus_musculus.GRCm38.102.gene.id_name.tsv", sep = "\t", header = T, check.names = F)
rownames(gene.table) = gene.table$`gene id`
gene.table.ambiguous = read.table("/media/biogenger/D/Reference/mouse/GRCm38_mm10/ensembl/diffID_sameNAME.tsv", sep = "\t", header = T, check.names = F)

count.m = as.matrix(readRDS(rds.mat)$umicount$inex$all)
barcode.df = read.csv(barcode.csv)
if (dim(barcode.df)[1] != dim(count.m)[2]) {
  warning(paste0("The barcode num of barcode_csv is ", dim(barcode.df)[1], ", not equal to the output matrix ", dim(count.m)[2], "!"))
}
rownames(barcode.df) = paste0(barcode.df$bc_B, barcode.df$bc_A)
barcode.df = barcode.df[colnames(count.m), ]
colnames(count.m) = paste0(barcode.df$iB, "x", barcode.df$iA)
print(paste("Input ID-matrix size is", dim(count.m)[1], "x", dim(count.m)[2]))

gene.table.ambiguous$hit = gene.table.ambiguous$`gene id` %in% rownames(count.m)
gene.table.ambiguous$expr = rowSums(count.m)[gene.table.ambiguous$`gene id`]
gene.table.ambiguous$expr[is.na(gene.table.ambiguous$expr)] = 0
for (gene in unique(gene.table.ambiguous$`gene name`)) {
  tmp = gene.table.ambiguous[gene.table.ambiguous$`gene name` %in% gene, ]
  tmp = tmp %>% filter(hit == TRUE)
  id.n = length(tmp$hit)
  hit.n = sum(tmp$hit)
  if (hit.n > 1) {
    print(tmp)
    print("Ids with the same gene name!")
    print(count.m[tmp$`gene id`, ][, 1:10])
    sum.expr = colSums(count.m[tmp$`gene id`, ])
    count.m[tmp$`gene id`, ] = 0
    count.m[tmp$`gene id`[1], ] = sum.expr
    print("Sum up ids with the same gene name to the first id and remove other ids!")
    print(count.m[tmp$`gene id`, ][, 1:10])
    print("---------------------------------------------------------------------------------------------------")
    count.m = count.m[!(rownames(count.m) %in% tmp$`gene id`[2:id.n]), ]
  }
}
rownames(count.m) = gene.table[rownames(count.m), "gene name"]

if (sum(duplicated(rownames(count.m))) > 0) {
  print("Duplicated gene names found!")
} else {
  print("Successfully!")
}
print(paste("Output Symbol-matrix size is", dim(count.m)[2], "x", dim(count.m)[1]))

write.table(t(count.m), file = name.mat, row.names = T, sep = "\t", quote = F)

