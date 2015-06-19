
args <- commandArgs(trailingOnly = TRUE)
f <- args[1]
load(f)
write.table(gwas.gerp, file=paste(f, ".txt", sep=""), quote=FALSE, row.names=FALSE)

