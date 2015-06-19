d <- read.table( "dat.txt", header=T)

top <- d[d$rs_number %in% c("rs11692564", "rs188303909", "rs6542457"),]
top <- top[order(top$rs_number),]
d <- d[!d$rs_number %in% c("rs11692564", "rs188303909", "rs6542457"),]
nonsig <- d[ d$p.value > 1.2e-5, ]
sig <- d[ d$p.value <= 1.2e-5, ]
print(length(sig$p.value))

postscript(file="panelb.ps",  width=3.5, height=1.75, pointsize=5, horizontal=FALSE, paper="special")
par(mar=c(4,5,2,0))
plot( sig$position, -log10(sig$p.value), 
      pch=21, col="black",bg="#0000ff", 
      ##cex=1.75, cex.axis=1.5, cex.lab=1.5,
      xaxt='n', 
      yaxt='n',
      #xlab="Position (bp/10^4)", 
      xlab="Position (\U00D7 10 kb)",
      ylab="LS BMD P value (-log10(P))", bty="n", ylim=c(0,-log10(1e-15)), 
      )

points(nonsig$position, -log10(nonsig$p.value), 
       pch=21, col="black",bg="#0000ff",
       ##cex=1.75, cex.axis=1.5, cex.lab=1.5
      )

abline(h=-log10(1.2e-8), lwd=1, lty=1, col="darkgreen")
points(top$position[seq(1,6,2)], -log10(top$p.value[seq(1,6,2)]), pch=17, col=rep("#ff0000",6))
points(top$position[seq(2,6,2)], -log10(top$p.value[seq(2,6,2)]), pch=19, col=rep("#ff0000",6))


axis(1, at=seq(119100000,119600000, by=100000), 
     labels=seq(119100000,119600000, by=100000)/10000, cex.axis=1, lwd=1)
axis(2, at=c(0,3,6,9,12,15), labels=c(0,3,6,9,12,15), cex.axis=1, las=2, lwd=1)
dev.off()