

args <- commandArgs(trailingOnly = TRUE)
gwas.f <- args[1]
gwas <- read.table(gwas.f, header=TRUE)
save(gwas, file=paste(gwas.f, ".Rdata", sep=""))
