################################################################################


library(skatMeta)

source("~/common/gefos.seq/SKATMETA/bin/skatMetaWrapper.devel.R")
args <- commandArgs(trailingOnly = TRUE)

region.f = args[1]
cohort.f = args[2]

source(cohort.f)

regions = read.table(region.f, header=T, as.is=T)

cwd = getwd()
setwd("~/common/gefos.seq/RESULTS/")
message("[**] RUNNING SKATMETA")
res <- skatMetaWrapper(cohorts, regions, skat.method="skatOMeta", skat.mafRange = list(c(0, 0.01)))
setwd(cwd)
message("[**] WRITING RESULTS")
write.table(res$meta, file="", quote=F, row.names=F)
message("[**] WRITING SNPINFO")
write.table(res$snpinfo,
            file=paste("out/", basename(region.f),
              ".snpinfo.windows.txt", sep=""),
            quote=F, row.names=F)
message("[**] DONE")
