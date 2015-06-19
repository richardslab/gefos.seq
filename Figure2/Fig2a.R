plot.it <- function(f, out.f, ph="Pheno", main="Title", y.extra=0.1, maf.cutoffs=c(0.005, 0.01, 0.05, 0.1, 0.5), power, boxfill, boxcol){
  
  d <- read.table(f, header=TRUE)
  d$beta <- abs(d$beta)
  
  ## bin 1
  d1 <- d[d$eaf<=0.01 & d$eaf>=0.005 & d$n_studies>2,]
  sb1 <- d1[which.min(d1$p.value),]$beta
  mb1 <- 0
  if (length(d1$p.value) > 0){
    mb1 <- mean(d1$beta)
  }else{
    sb1 <- NA
  }
  ## bin 2
  d2 <- d[d$eaf<=0.05 & d$eaf>0.01 & d$n_studies>2,]
  mb2 <- mean(d2$beta)
  sb2 <- d2[which.min(d2$p.value),]$beta
  ## bin 3
  d3 <- d[d$eaf<=0.10 & d$eaf>0.05 & d$n_studies>2,]
  mb3 <- mean(d3$beta)
  sb3 <- d3[which.min(d3$p.value),]$beta
  ## bin 4
  d4 <- d[d$eaf>0.10 & d$n_studies>2,]
  mb4 <- mean(d4$beta)
  sb4 <- d4[which.min(d4$p.value),]$beta
  
  par(bty="n", xaxt="n")
  n <- power[power$ph == ph,]$n[1]
  boxplot(d1$beta,d2$beta,d3$beta,d4$beta, col=boxfill, border=boxcol, bty="n",
          pars = list(bty="n"), cex=1, cex.main=1, lwd=0.8,
          outline=TRUE, pch=20, ylab="Absolute Value of Beta", main=main, 
          ylim=c(0, max(c(d1$beta,d2$beta,d3$beta,d4$beta)+y.extra))
          #ylim=c(0,0.5)
  )

  i <- 1
  for (betas in list(d1$beta, d2$beta, d3$beta, d4$beta)){
    barh <- min(c(betas, power[power$maf==maf.cutoffs[i+1] & power$ph == ph,]$beta))
    rect(i-0.4, 
         0,#power[power$maf==maf.cutoffs[i] & power$ph == ph,]$beta,  
         i+0.4, 
         barh,
         col="#cccccc", border="#ffffff", lwd=0)
    i <- i + 1
  }
  
  title(xlab="MAF Range (%)", line=2)
  par(xaxt="s")
  axis(1, at=c(1,2,3,4), labels=c("0.5-\n1", "1-\n5", "5-\n10", "10-\n50"), tick=FALSE, 
       line=-1, cex.axis=0.75)        
  
}

power = read.table("power.txt", header=TRUE)
maf.cutoffs <- c(0.005, 0.01, 0.05, 0.1, 0.5)

postscript(file="boxplots.1.2e-8.2.ps", width=3.75, height=1.5, paper="special", horizontal = FALSE, pointsize=5)
par(mfrow=c(1,3), cex.axis=1, cex.main=1, cex.lab=1, las=1, mar=c(3,4,2,2), ask=F)
par(bty="n")
plot.it("fa2stu.MAF0.005.1.2e-8.txt", "fa2stu.out.jpg", "FA", "Forearm", 0.1, maf.cutoffs, power, "#ffcccc", "#ff0000")
plot.it("fn2stu.MAF0.005.1.2e-8.txt", "fn2stu.out.jpg", "FN", "Femoral Neck", 0.2, maf.cutoffs, power, "#006600", "#008800")
plot.it("ls2stu.MAF0.005.1.2e-8.txt", "ls2stu.out.jpg", "LS", "Lumbar Spine", 0.07, maf.cutoffs, power, "#0000ff", "#000033")
dev.off()

postscript(file="boxplots.pruned.1.2e-8.2.ps",  width=3.75, height=1.5, paper="special", horizontal = FALSE, pointsize=7)
par(mfrow=c(1,3), las=1, mar=c(3,4,2,2), ask=F, oma=c(0,5,0,5))
par(bty="n")
plot.it("fa2stu.out.pruned.txt", "fa2stu.out.jpg", "FA", "Forearm",
        0.1, maf.cutoffs, power, "#ffcccc", "#ff0000")
plot.it("fn2stu.out.pruned.txt", "fn2stu.out.jpg", "FN", "Femoral Neck",
        0.2, maf.cutoffs, power, "#ccffcc", "#006600")
plot.it("ls2stu.out.pruned.txt", "ls2stu.out.jpg", "LS", "Lumbar Spine",
        0.07, maf.cutoffs, power, "#ccccff", "#0000ff")
dev.off()

