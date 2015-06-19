for (i in seq(0,4,by=0.2)) {
  f <- paste(i,".Rdata", sep="")
  if (file.exists(f)){
    load(f)
    sm.out$sum.all$bin <- i
    sm.out$sum.matched$bin <- i
    head=F
    if (i == 0){
      head=T
    }
    write.table(rbind(sm.out$sum.all, sm.out$sum.matched), file="", sep=",", quote=FALSE, col.names=head)
  }
}
