
write("file treat.prop treat.n treat.se treat.below.q treat.fdr.n contr.prop contr.n contr.se cont.below.q cont.fdr.n", file="")
for (i in seq(0,4,by=0.2)) {
  f <- paste(i,".Rdata", sep="")
  if (file.exists(f)){
    load(f)
    treat.n <- treatment.fdr$n
    treat.se <- sqrt((treatment.fdr$prop * (1 - treatment.fdr$prop)) / treatment.fdr$n)
    contr.n <- control.fdr$n
    contr.se <- sqrt((control.fdr$prop * (1 - control.fdr$prop)) / control.fdr$n)
    write(paste(f, treatment.fdr$prop, treat.n, treat.se, treatment.fdr$below.q,  treatment.fdr$n, control.fdr$prop, contr.n, contr.se, control.fdr$below.q, control.fdr$n, sep=" "), file="")
  }
}
