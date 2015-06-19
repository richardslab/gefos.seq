
d <- read.csv("data/microct.csv")

postscript(file="microct.3.fig3.new.ps", width=3.5, height=2, pointsize=5, paper="special", horizontal = FALSE)
par(mar=c(2,4,0,1), cex=1.3)

# SPINE L-5
spine <- d[(d$site == "spine") & (d$vert=="L-5"),]
cort <- d[(d$site == "fem") & (d$vert=="cort"),]

x1 <- 1
x2 <- 2

g1 <- spine[(spine$genotype=="En1lox/+"),]$Tb.N
g2 <- spine[(spine$genotype=="En1Cre/lox"),]$Tb.N
total <- sum(c(g1, g2))
g1 <- g1/rep(total, length(g1))*100
g2 <- g2/rep(total, length(g2))*100
r1 <- t.test(g2, g1, paired = FALSE)
g12 <- c(g1,g2)
ymax <- ceiling(max(g12))
ymin <- floor(min(g12))

ylab=""
plot(0,0,xlim=c(0,12), ylim=c(4,18), 
     xaxt="n", yaxt="n", xlab="", 
     col=0, main="", font.main=1, bty="n", ylab="% of Total")
axis(2, at=seq(4,18,4),
     labels=seq(4,18,4), las=2)

kocol = "#660000"
wtcol = "#006600"
kobg = "#cc0000"
wtbg = "#00cc00"
points(c(rep(x1,length(g1))), c(g1), col=wtcol, pch=21, lwd=1, bg=wtbg)
points(c(rep(x2,length(g2))), c(g2), col=kocol, pch=21, lwd=1, bg=kobg)
lines(c(x1-0.4, x1+0.4), rep(mean(g1, na.rm=TRUE),2), col=wtcol, lwd=1)
lines(c(x2-0.4, x2+0.4), rep(mean(g2, na.rm=TRUE),2), col=kocol, lwd=1)

x1 <- 4
x2 <- 5
g1 <- spine[(spine$genotype=="En1lox/+"),]$BV.TV
g2 <- spine[(spine$genotype=="En1Cre/lox"),]$BV.TV
total <- sum(c(g1, g2))
g1 <- g1/rep(total, length(g1))*100
g2 <- g2/rep(total, length(g2))*100
r2 <- t.test(g2, g1, paired = FALSE)
g12 <- c(g1,g2)
ymax <- ceiling(max(g12))
ymin <- floor(min(g12))
points(c(rep(x1,length(g1))), c(g1), col=wtcol, pch=21, lwd=1, bg=wtbg)
points(c(rep(x2,length(g2))), c(g2), col=kocol, pch=21, lwd=1, bg=kobg)
lines(c(x1-0.4, x1+0.4), rep(mean(g1, na.rm=TRUE),2), col=wtcol, lwd=2)
lines(c(x2-0.4, x2+0.4), rep(mean(g2, na.rm=TRUE),2), col=kocol, lwd=2)

x1 <- 7
x2 <- 8
g1 <- spine[(spine$genotype=="En1lox/+"),]$Tb.Th
g2 <- spine[(spine$genotype=="En1Cre/lox"),]$Tb.Th
total <- sum(c(g1, g2))
g1 <- g1/rep(total, length(g1))*100
g2 <- g2/rep(total, length(g2))*100
r3 <- t.test(g2, g1, paired = FALSE)
g12 <- c(g1,g2)
ymax <- ceiling(max(g12))
ymin <- floor(min(g12))
points(c(rep(x1,length(g1))), c(g1), col=wtcol, pch=21, lwd=1, bg=wtbg)
points(c(rep(x2,length(g2))), c(g2), col=kocol, pch=21, lwd=1, bg=kobg)
lines(c(x1-0.4, x1+0.4), rep(mean(g1, na.rm=TRUE),2), col=wtcol, lwd=1)
lines(c(x2-0.4, x2+0.4), rep(mean(g2, na.rm=TRUE),2), col=kocol, lwd=1)


d2 <- read.csv("data/Joyner2V_summaryOnly.csv", header=TRUE)
x1 <- 10
x2 <- 11
g1 <- d2[(d2$Bone.type=="Control"),]$TRAP_BS
g2 <- d2[(d2$Bone.type=="Test"),]$TRAP_BS
total <- sum(c(g1, g2))
g1 <- g1/rep(total, length(g1))*100
g2 <- g2/rep(total, length(g2))*100
r4 <- t.test(g2, g1, paired = FALSE)
g12 <- c(g1,g2)
ymax <- ceiling(max(g12))
ymin <- floor(min(g12))
points(c(rep(x1,length(g1))), c(g1), col=wtcol, pch=21, lwd=1, bg=wtbg)
points(c(rep(x2,length(g2))), c(g2), col=kocol, pch=21, lwd=1, bg=kobg)
lines(c(x1-0.4, x1+0.4), rep(mean(g1, na.rm=TRUE),2), col=wtcol, lwd=1)
lines(c(x2-0.4, x2+0.4), rep(mean(g2, na.rm=TRUE),2), col=kocol, lwd=1)

axis(1, at=c(1.5,4.5,7.5,10.5),
     labels=paste(c("Tb.N","BV/TV","Tb.Th", "TRAP/BS"),  sprintf("%0.1e", c(r1$p.value, r2$p.value, r3$p.value, r4$p.value)), sep="\n"),
     tick=FALSE,
     font=1,
     cex.axis=1)

legend(-0.5,18,c("Control", "sdEn1"), bty="n", pch=c(19,19), col=c(wtbg, kobg), cex=0.6)
dev.off()
