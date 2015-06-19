ccol = "black"
tcol = "black"
cpch = 21
tpch = 19

a <- read.table("table.txt", header=T)
x <- 1:3
#jpeg("csq.jpg", quality=100)
par(mar=c(12,12,6,12))
plot(x, a$csq.prop, xlim=c(0.7,3), ylim=c(0.0000,0.0004), col=ccol, pch=cpch, main="Coding", xlab="", ylab="Proportion below FDR q-value of 0.05", bty="n", xaxt="n")
points(x, a$cds.prop, col=tcol, pch=tpch)
axis(1, at=1:3, labels=c("Forearm", "Lumbar Spine", "Femoral Neck"), las=2)

for (i in 1:length(x)){
  print(x[i])
  lines(c(x[i], x[i]), c( a[i,]$csq.prop - a[i,]$csq.se,  a[i,]$csq.prop + a[i,]$csq.se), col=ccol)
  lines(c(x[i], x[i]), c( a[i,]$cds.prop - a[i,]$cds.se,  a[i,]$cds.prop + a[i,]$cds.se), col=tcol)
}
#
#dev.off()
