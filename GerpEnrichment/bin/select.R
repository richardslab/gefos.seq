library(plyr)
if (TRUE){
message("LOADING META")
load("../../gene_distance/fdr.fa/fa2stu.gerp.gene.bed.trim.space.meta.trim.Rdata")
colnames(d) <- c("chrom", "start", "end", "gerp", "dist", "rs", "eaf", "pv")
message("LOADING R2")
r2.file <- "22.fa2stu.maf0.005.ld.tab"
r2 <- read.table(r2.file, header=T)

message("MERGING")
r2.d <- merge(r2, d, by.x="SNP_B", by.y="rs", all.x=TRUE)
rm(d)
r2.d2 <- r2.d[,c("SNP_A", "SNP_B", "pv")]
}
message("PICKING MINP")
r2.d2.minp <- ddply(r2.d2, .(SNP_A), subset, pv==min(pv))
