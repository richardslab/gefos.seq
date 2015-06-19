
args <- commandArgs(trailingOnly = TRUE)

gwas.f <- args[1]
gerp.f <- args[2]
chrom <- args[3]

message("Loading")
load(gwas.f)
load(gerp.f)

message("Aggregating GERP")
gerp <- aggregate(gerp$scores, list(gerp$bp), mean)
colnames(gerp) <- c("bp", "scores")

message("Merging")
gwas.gerp <- merge(gwas, gerp, by.x="position", by.y="bp", all.x=TRUE, all.y=FALSE)

message("Saving")
save(gwas.gerp, file=paste(chrom, ".meta.gerp.Rdata", sep=""))
