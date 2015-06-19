load.RDS <- function (f) { return(readRDS(file=f)) }
load.df <- function(f) { return(read.table(f, header=T)) }

cohort.phenotype = "LS"

cohorts <- list(
             "TWINSUK.seq" = list(name="TWINKSUK.seq", type="genome", samples="ALL.A",
               skatcohort.file=function (chrom, gene) {
                 return(paste("TWINSUK/genome-seq.skatcohort/", cohort.phenotype, "/ALL.A/gefos/",
                              chrom, "_", gene, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom, gene){
                 return(paste("TWINSUK/genome-seq.skatcohort/snpinfo/",
                              chrom, "_", gene, ".snpinfo__clean__", sep=""))},
               snpinfo.load=load.df
               ),
             "TWINSUK.imp" = list(name="TWINKSUK.imp", type="genome", samples="ALL.A",
               skatcohort.file=function (chrom, gene) {
                 return(paste("TWINSUK/genome-imp.skatcohort/", cohort.phenotype, "/ALL.A/gefos/",
                              chrom, "_", gene, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom, gene){
                 return(paste("TWINSUK/genome-imp.skatcohort/snpinfo/",
                              chrom, "_", gene, ".snpinfo__clean__.2ref", sep=""))},
               snpinfo.load=load.df
               ),
             "ALSPAC.imp" = list(name="ALSPAC.imp", type="genome", samples="ALL.A",
               skatcohort.file=function (chrom, gene) {
                 return(paste("ALSPAC/genome-imp.skatcohort/", cohort.phenotype, "/ALL.A/",
                              chrom, "_", gene, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom, gene){
                 return(paste("ALSPAC/genome-imp.skatcohort/snpinfo/",
                              chrom, "_", gene, ".snpinfo__clean__.2ref", sep=""))},
               snpinfo.load=load.df
               ),
             "ALSPAC.seq" = list(name="ALSPAC.seq", type="genome", samples="ALL.A",
               skatcohort.file=function (chrom, gene) {
                 return(paste("ALSPAC/genome-seq.skatcohort/", cohort.phenotype, "/ALL.A/",
                              chrom, "_", gene, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom, gene){
                 return(paste("ALSPAC/genome-seq.skatcohort/snpinfo/",
                              chrom, "_", gene, ".snpinfo__clean__", sep=""))},
               snpinfo.load=load.df
               ),
             "AOGC.imp" = list(name="AOGC.imp", type="genome", samples="FEMALES.A-b",
               skatcohort.file=function (chrom, gene) {
                 return(paste("AOGC/genome-imp.skatcohort/", cohort.phenotype, "/FEMALES.A-b/gefos/",
                              chrom, "/", chrom, "_", gene, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom, gene){
                 return(paste("AOGC/genome-imp.skatcohort/", cohort.phenotype, "/FEMALES.A-b/snpinfo/chr",
                              chrom, "/snpinfo/", chrom, "_", gene, ".snpinfo__clean__.2ref",
                              sep=""))},
               snpinfo.load=load.df
               ),
             "FHAM.imp" = list(name="FHAM.imp", type="genome", samples="ALL.A-b",
               skatcohort.file=function (chrom, gene) {
                 return(paste("FRAMINGHAM/genome-imp.skatcohort/", cohort.phenotype, "/ALL.A-b/",
                              chrom, "_", gene, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom, gene){
                 return(paste("FRAMINGHAM/genome-imp.skatcohort/snpinfo/",
                              chrom, "_", gene, ".snpinfo__clean__.2ref", sep=""))},
               snpinfo.load=load.df
               ),
             "FHAM.exo" = list(name="FHAM.exo", type="exome", samples="ALL.C",
               skatcohort.file=function (chrom) {
                 return(paste("FRAMINGHAM/exome-seq.skatcohort/", cohort.phenotype, "/ALL.C/",
                              chrom, ".rds__fixed__", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom){
                 return(paste("FRAMINGHAM/exome-seq.skatcohort/", cohort.phenotype, "/ALL.C/",
                              chrom, ".snpinfo__fixed____poly__", sep=""))},
               snpinfo.load=load.df
               ),
             "AOGC.exo" = list(name="AOGC.exo", type="exome", samples="FEMALES.C",
               skatcohort.file=function (chrom) {
                 return(paste("AOGC/exome-seq.skatcohort/", cohort.phenotype, "/FEMALES.C/",
                              chrom, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom){
                 return(paste("AOGC/exome-seq.skatcohort/", cohort.phenotype, "/FEMALES.C/",
                              chrom, ".snpinfo__poly__", sep=""))},
               snpinfo.load=load.df
               ),
             "WHI.exo" = list(name="WHI.exo", type="exome", samples="ALL.C",
               skatcohort.file=function (chrom) {
                 return(paste("WHI/exome-seq.skatcohort/", cohort.phenotype, "/ALL.C/",
                              chrom, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom){
                 return(paste("WHI/exome-seq.skatcohort/", cohort.phenotype, "/ALL.C/",
                              chrom, ".snpinfo__poly__", sep=""))},
               snpinfo.load=load.df
               ),
             "MROS.imp" = list(name="MROS.imp", type="genome", samples="ALL.A",
               skatcohort.file=function (chrom, gene) {
                 return(paste("MROS/genome-imp.skatcohort/", cohort.phenotype, "/ALL.A/",
                              chrom, "_", gene, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom, gene){
                 return(paste("MROS/genome-imp.skatcohort/snpinfo/",
                              chrom, "_", gene, ".snpinfo__clean__.2ref", sep=""))},
               snpinfo.load=load.df
               ),
             "SOF.imp" = list(name="SOF.imp", type="genome", samples="ALL.A",
               skatcohort.file=function (chrom, gene) {
                 return(paste("SOF/genome-imp.skatcohort/", cohort.phenotype, "/ALL.A/",
                              chrom, "_", gene, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom, gene){
                 return(paste("SOF/genome-imp.skatcohort/snpinfo/",
                              chrom, "_", gene, ".snpinfo__clean__.2ref", sep=""))},
               snpinfo.load=load.df
               ),
             "ROTTERDAM.imp" = list(name="ROTTERDAM.imp", type="genome", samples="ALL.A-b",
               skatcohort.file=function (chrom, gene) {
                 return(paste("ROTTERDAM/genome-imp.skatcohort/", cohort.phenotype, "/ALL.A-b/",
                              chrom, "_", gene, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom, gene){
                 return(paste("ROTTERDAM/genome-imp.skatcohort/snpinfo/",
                              chrom, "_", gene, ".snpinfo__clean__.2ref", sep=""))},
               snpinfo.load=load.df
               ),
             "ROTTERDAM.exo" = list(name="ROTTERDAM.exo", type="exome", samples="ALL.C",
               skatcohort.file=function (chrom) {
                 return(paste("ROTTERDAM/exome-seq.skatcohort/", cohort.phenotype, "/ALL.C/",
                              chrom, ".rds", sep=""))},
               skatcohort.load=load.RDS,
               snpinfo.file=function (chrom){
                 return(paste("ROTTERDAM/exome-seq.skatcohort/", cohort.phenotype, "/ALL.C/",
                              chrom, ".snpinfo__poly__", sep=""))},
               snpinfo.load=load.df
               )
             )
