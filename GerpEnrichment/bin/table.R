
write("file treat.prop treat.n treat.se treat.below.q treat.fdr.n contr.prop contr.n contr.se cont.below.q cont.fdr.n", file="")
for (i in seq(0,4,by=0.2)) {
  f <- paste(i,".Rdata", sep="")
  if (file.exists(f)){
    load(f)
    treat.n <- sm.out$nn[[6]]
    treat.se <- sqrt((treatment.fdr$prop * (1 - treatment.fdr$prop)) / treat.n)
    contr.n <- sm.out$nn[[2]]
    contr.se <- sqrt((control.fdr$prop * (1 - control.fdr$prop)) / contr.n)
    write(paste(f, treatment.fdr$prop, treat.n, treat.se, treatment.fdr$below.q,  treatment.fdr$n, control.fdr$prop, contr.n, contr.se, control.fdr$below.q, control.fdr$n, sep=" "), file="")
  }
}
