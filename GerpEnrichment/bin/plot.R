


plot.it <- function(ph, main) {
  
  cpch = 21
  tpch = 19
  table.f = paste(ph, "/table.txt", sep="")
  a <- read.table(table.f, header=T)
  a <- a[2:(length(a[,1])-1),]
  x <- seq(0.2,3.8,by=0.2)
  plot(x, a$treat.prop, ylim=c(0, 0.0035), col=ccol, pch=cpch, main=main, xlab="GERP threshold", ylab="", bty="n", xaxt="n")
  points(x, a$contr.prop, col=tcol, pch=tpch)
  axis(1, at=x)
  
  for (i in 1:length(x)){
    print(x[i])
    lines(c(x[i], x[i]), c( a[i,]$treat.prop - a[i,]$treat.se,  a[i,]$treat.prop + a[i,]$treat.se), col=ccol)
    lines(c(x[i], x[i]), c( a[i,]$contr.prop - a[i,]$contr.se,  a[i,]$contr.prop + a[i,]$contr.se), col=tcol)
  }
}

phenos <- list(list("ls", "Lumbar Spine"),
          list("fn", "Femoral Neck"),
          list("fa", "Forearm"))

png(file="figure.png", res=240, width=2400, height=1200)


tcol = "#ff3333"
ccol = "#006600"
cpch = 21
tpch = 19
par(mfrow=c(1,5), mar=c(8,4,4,2))
a <- read.table("table.syn.txt", header=T)
x <- 1:3
plot(x, a$csq.prop, xlim=c(0.7,3), ylim=c(0.0000,0.0004), col=ccol, pch=cpch, main="Synonymous", xlab="", ylab="Proportion below FDR q-value of 0.05", bty="n", xaxt="n")
points(x, a$cds.prop, col=tcol, pch=tpch)
axis(1, at=1:3, labels=c("Forearm", "Lumbar Spine", "Femoral Neck"), las=2)

for (i in 1:length(x)){
  print(x[i])
  lines(c(x[i], x[i]), c( a[i,]$csq.prop - a[i,]$csq.se,  a[i,]$csq.prop + a[i,]$csq.se), col=ccol)
  lines(c(x[i], x[i]), c( a[i,]$cds.prop - a[i,]$cds.se,  a[i,]$cds.prop + a[i,]$cds.se), col=tcol)
}
par(mar=c(8,2,4,2))
cpch = 21
tpch = 19

a <- read.table("table.csq.txt", header=T)
x <- 1:3
plot(x, a$csq.prop, xlim=c(0.7,3), ylim=c(0.0000,0.0004), col=ccol, pch=cpch, main="Deleterious", xlab="", ylab="", bty="n", xaxt="n")
points(x, a$cds.prop, col=tcol, pch=tpch)
axis(1, at=1:3, labels=c("Forearm", "Lumbar Spine", "Femoral Neck"), las=2)

for (i in 1:length(x)){
  print(x[i])
  lines(c(x[i], x[i]), c( a[i,]$csq.prop - a[i,]$csq.se,  a[i,]$csq.prop + a[i,]$csq.se), col=ccol)
  lines(c(x[i], x[i]), c( a[i,]$cds.prop - a[i,]$cds.se,  a[i,]$cds.prop + a[i,]$cds.se), col=tcol)
}


for (ph in phenos){
  plot.it(ph[[1]][[1]], ph[[2]][[1]])
}

dev.off()
