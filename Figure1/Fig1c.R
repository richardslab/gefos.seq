d <- read.table("snps.tbl", header=TRUE, sep="\t")

postscript(file="panelc.ps", width=1.75, height=1.75, pointsize=5, horizontal=FALSE, paper="special")
par(mar=c(4,5,2,2))

plot(d$Freq1, abs(d$LS.S2.Beta), xlim=c(0,1), ylim=c(0,0.25), pch=21, bty="l",
     xlab="Allele Frequency of Genome-wide\nSignificant Variants",
     ylab="Absolute Effect Size on BMD\nin Standard Deviations",
     las=1, col="black", bg="blue",
     lwd=0.5)
abline(h=mean(abs(d$LS.S2.Beta)), col="#ff0000", lty=1, lwd=1, xpd=FALSE)
gefos.maf = c(0.017, 0.022, 0.053)
gefos.effect = c(0.20, 0.14, 0.09)

points(gefos.maf, gefos.effect, col="black", pch=21, bg="red", lwd=0.8)

dev.off()
