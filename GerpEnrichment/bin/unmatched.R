library(MatchIt)
library(fdrtool)

fdrtool.cutoff <- function(pv, alpha=0.05) {
  obj <- fdrtool(pv, statistic="pvalue", plot=F, verbose=FALSE)
  qv <- obj$qval
  pi0 <- obj$param[1, "eta0"]
  indx <- (qv < alpha)
  p1 <- rep(0, length(pv))  
  p1[indx] <- pv[indx]
  p2 <- rep(1, length(pv))
  p2[!indx] <- pv[!indx]
  pvcutoff <- (max(p1) + min(p2))/2
  below.q <- length(which(qv < alpha))
  n <- length(qv)
  prop <- length(which(qv < alpha))/length(qv)
  r <- list(p.cutoff=pvcutoff, prop=prop, below.q=below.q, n=n)
  r
}

args <- commandArgs(trailingOnly=TRUE)
gerp.threshold <- as.numeric(args[1])
meta.file <- args[2]

#gerp.threshold <- 1.2
message("Loading data ...")
## d <- read.table("ls.100k.meta", header=F, sep=" ")
## save(d, file="ls2stu.gerp.gene.bed.trim.space.meta.trim.Rdata")
load(meta.file)
colnames(d) <- c("chrom", "start", "end", "gerp", "dist", "rs", "eaf", "pv")
message("Cleaning ...")
d <- d[which(!is.na(d$gerp)),]
d$treat <- as.integer(d$gerp >= gerp.threshold)

Y <- d$pv
Tr <- d$treat
X <- cbind(d$dist, d$eaf)

message("Matching")
#m.out <- matchit(Tr ~ eaf + dist, data = d, method="cem")
#m.d <- match.data(m.out)

treatment.fdr <- fdrtool.cutoff(d$pv[ d$treat == 1 ])
message("FDR Controls")
control.fdr <- fdrtool.cutoff(d$pv[ d$treat == 0 ])

out.file = paste(gerp.threshold, ".Rdata", sep="", collapse="")

#sm.out <- summary(m.out)
save(gerp.threshold, treatment.fdr, control.fdr, file=out.file)
