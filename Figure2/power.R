power<- function(beta, sigma=1, R2=1, p, n, alpha)
{
  ## simple linear regression: (\hat\beta)^2/Var(\hat\beta) ~ chisq(df=1, ncp)
  ## Signal-to-noise ratio: snr = (beta/sigma)^2
  snr<- (beta/sigma)^2
  Qa<- qchisq(1 - alpha, df=1, ncp=0)
  NCP<- 2*(n-1)*p*(1-p)*snr*R2
  pw<- pchisq(Qa, df=1, ncp=NCP, lower.tail = FALSE)
  names(pw)<- alpha
  pw
}

cohorts <- list(
             list(ph = "FA", n = 7822),
             list(ph = "LS", n = 25195),
             list(ph = "FN", n = 29156)
             )

alpha <- 5e-8
maf.list <- c(0.005, 0.01, 0.05, 0.10, 0.5)
#maf.list <- c(0.01)
out.f <- "power.txt"
write("beta maf n alpha ph power", file=out.f)
min.power <- 0.8
for (cohort in cohorts){
  for (p in maf.list){
    for (b in seq(0,1,0.001)){
      pw <- power(beta=b, p=p, n=cohort$n, alpha=5e-8)
      if (pw>=min.power) {
        write(paste(b, p, cohort$n, alpha=alpha, cohort$ph, as.numeric(pw), sep=" "), file=out.f, append=TRUE)
        break;
      }
    }
  }
}

