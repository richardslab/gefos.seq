#!/usr/bin/R
#################################################################################
## CAVEATS:
##    1. Drop-one analysis not supported for singlesnpMeta

library(skatMeta)
options(warn = 1)
## append item to list
lappend <- function(lst, obj) {
  lst[[length(lst) + 1]] <- obj
  return(lst)
}

remove.single.cohort.singletons <- function(snpinfo.list, skatcohort.list){
  duplicates <- which(duplicated(snpinfo.list))
  
  single.cohort.snps <- snpinfo.list[!snpinfo.list$Name %in%
                                     snpinfo.list$Name[duplicates],]$Name
  duplicated.cohort.snps <- snpinfo.list[snpinfo.list$Name %in%
                                         snpinfo.list$Name[duplicates],]$Name
  message("Total cohort SNP count = ", length(snpinfo.list$Name))
  message("Single cohort SNP count = ", length(single.cohort.snps))
  message("Multiple cohort SNP count = ", length(duplicated.cohort.snps))
  
  singletons = c()
  for (snp in single.cohort.snps){
    for (skatcohort.item in skatcohort.list){
      cohort.snps = colnames(skatcohort.item[[1]]$cov)
      cohort.name = skatcohort.item[[1]]$cohortName
      if (snp %in% cohort.snps){
        maf = skatcohort.item[[1]]$maf[which(snp == cohort.snps)]
        alt.count = maf * skatcohort.item[[1]]$n
        #message(maf, " ", alt.count, " ", snp, " ", cohort.name)
        if (alt.count < 2) {
          singletons <- c(singletons, snp)
        }
      }
    }
    
    
  }
  message("Single Cohort Singletons = ", length(singletons))
  single.cohort.snps = single.cohort.snps[which(!single.cohort.snps %in% singletons)]
  message("Single cohort SNP count, singletons removed = ", length(single.cohort.snps))
  cohort.snps = c(single.cohort.snps, duplicated.cohort.snps)
  snpinfo.list.clean <- snpinfo.list[snpinfo.list$Name %in% cohort.snps,]
  message("Total SNPs retained = ", length(snpinfo.list.clean$Name))
  return(snpinfo.list.clean)
}

get.positional.windows <- function(snpinfo.list, gene.start, gene.end,
                                   positional.window.size,
                                   positional.window.step, gene, log.tag) {
  snpinfo.window.list <- list()
  win.id <- 1
  for (win.s in seq(gene.start, gene.end, by=positional.window.step)){
    win.e <- win.s + positional.window.size - 1
    if (win.e > gene.end) {
      win.e <- gene.end
    }
    region <- paste(gene, ".", win.s, "-", win.e, sep="")
    snpinfo.window <- list(region = region, window.start = win.s,
                           window.end = win.e, window.id = win.id,
                           snpinfo = snpinfo.list[(snpinfo.list$Name >= win.s
                             & snpinfo.list$Name < win.e),])
    
    if (length(snpinfo.window$snpinfo$Name) == 0) {
      warning(log.tag, "WINDOW ID = ", win.id," HAS 0 SNPS")
      win.id <- win.id + 1
      next;
    }
    win.id <- win.id + 1
    snpinfo.window$snpinfo$window <- region
    snpinfo.window.list <- lappend(snpinfo.window.list, snpinfo.window)
    
  }
  return(snpinfo.window.list)
}

get.snp.windows <- function(snpinfo.list, window.size, window.step, nsnps, gene) {
  snpinfo.window.list <- list()
  win.id <- 1
  for (win.s in seq(1, nsnps, by=window.step)){
    skatmeta.params <- list()
    win.e <- win.s + window.size - 1
    if (win.e > nsnps){
      win.e <- nsnps
    }
    if ((win.e - win.s + 1) < window.size && (win.e > window.size)){
      win.s <- win.e - window.size + 1
    }
    window.start <- snpinfo.list[win.s,]$Name
    window.end <- snpinfo.list[win.e,]$Name
    region <- paste(gene, ".", win.s, "-", win.e, sep="")
    snpinfo.window <- list(region = region, window.start = window.start,
                           window.end = window.end, window.id = win.id,
                           snpinfo = snpinfo.list[win.s:win.e,])
    snpinfo.window$snpinfo$window <- region
    snpinfo.window.list <- lappend(snpinfo.window.list, snpinfo.window)
    win.id <- win.id + 1
    if (win.e == nsnps){
      return(snpinfo.window.list)
    }
  }
  return(snpinfo.window.list)
}

##' Run skatMeta analysis for GEFOS consortium data.
##'
##' Prepares data for meta-analysis including:
##'   * Loading skatCohort objects either by entire chromosome or gene.
##'   * Manipulate skatCohort objects: fetch gene data from chromosome object,
##'     filter SNPs from individual cohorts, etc.
##'   * Run meta-analysis, splitting regions into windows, and optionally
##'     drop-one meta-analysis for SNPs and cohorts
##' 
##' @title Run skatMeta analysis for GEFOS consortium data
##' @param cohorts An object containing cohorts to meta-analyze
##' @param regions An object containing regions to meta-analyze
##' @param window.size SNP window size
##' @param window.step SNP window step size
##' @param skat.method Method for SKAT (skatMeta, skatOMeta)
##' @param skat.mafRange MAF range for SKAT
##' @param skat.pvalue.method P-value method for SKAT
##' @param skat.meta.wts Weight function for SKAT.
##' @param do.dropone Boolean flag to enable/disable drop-one analysis
##' @param dropone.max.pvalue If drop-one enabled, threshold to perform drop-one
##' analysis
##' @return Object of meta-analysis results and object of SNPInfo for all windows
##' @author Vince Forgetta

skatMetaWrapper <- function(cohorts = cohorts, regions = regions,
                            window.size = 30, window.step = 20,
                            skat.method = "skatOMeta",
                            skat.mafRange = list(c(0, 0.05)),
                            skat.pvalue.method = "liu",
                            skat.meta.wts = function(maf) { dbeta(maf, 1, 25) },
                            do.dropone = FALSE,
                            dropone.max.pvalue = 5e-6,
                            do.positional.windows = FALSE,
                            positional.window.size = 1600,
                            positional.window.step = 1000,
                            do.remove.single.cohort.singletons = FALSE
                            ){
  ## Init results data.frame
  meta.res <- data.frame(
                region=character(), gene=character(),
                parent_region=character(), 
                #skat
                p.value=numeric(), qmeta=numeric(),
                cmaf=numeric(), nmiss=numeric(), nsnps=numeric(),
                #singlesnp
                maf=numeric(), name=character(), ntotal=numeric(),
                #burden
                beta=numeric(), se=numeric(),  cmafTotal=numeric(),
                cmafUsed=numeric(), nsnpsTotal=numeric(), nsnpsUsed=numeric(),
                #skatO
                rho=numeric(), pmin=numeric(), 
                errflag=numeric(),
                #custom
                nsubjects=character(), chrom=numeric(), window.start=numeric(),
                window.end=numeric(), widnow.id=numeric(),
                window.span=numeric(), window.nsnps=numeric(),  
                dropped.item=character(), cohorts=character(),
                stringsAsFactors=FALSE)
  
  ##
  ## PROCESS REGIONS BY CHROMOSOME
  ##

  snpinfo.all.windows <- data.frame(Name=character(), gene=character(),
    window=character())
  log.tag <- "[P] "
  message(log.tag, "Window.Size = ", window.size)
  message(log.tag, "Window.Step = ", window.step)
  message(log.tag, "Skat.Method = ", skat.method)
  message(log.tag, "Skat.MafRange = ", skat.mafRange)
  message(log.tag, "Skat.Pvalue.Method = ", skat.pvalue.method)
  message(log.tag, "Skat.Meta.Wts = ", gsub("\n", " ", list(skat.meta.wts)))
  message(log.tag, "Do.Dropone = ", do.dropone)
  message(log.tag, "Dropone.Max.Pvalue = ", dropone.max.pvalue)
  
  for (chrom in 1:22){
    
    ### Init data structs to collect info for chromosome-level skatMeta results
    skatcohort.exome.list <- list()
    snpinfo.exome.list <- list()
    log.tag <- paste("[",chrom,"] ", sep="")
    message(log.tag, "Processing Chromosome")
    ### Get regions for iteration's chromosome
    regions.chr <- regions[regions$CHR == chrom,]
    num.regions <- length(regions.chr[,1])
    message(log.tag, "Regions On Chromosome = ", num.regions)

    if (num.regions == 0){
      message(log.tag, "No Regions Found For Chromosome ... Skipping")
      next;
    }

    ##
    ### LOAD CHROMOSOME-LEVEL DATA, TAGGED AS "EXOME"
    ##
    exome.cohorts.loaded <- c()
    for (cohort in cohorts){
      ## If cohort data is provided by chromosome ("exome")
      if (cohort$type == "exome"){
        ## Check files exist
        if (file.exists(cohort$skatcohort.file(chrom)) &&
            file.exists(cohort$snpinfo.file(chrom))) {
          message(log.tag, "Loading Cohort ", cohort$name)
          ## Load skatCohort and SNPInfo data
          sk <- try(cohort$skatcohort.load(cohort$skatcohort.file(chrom)))
          si <- try(list(name=cohort$name,
                         snpinfo=cohort$snpinfo.load(cohort$snpinfo.file(chrom))))
          ## Check loading worked, append to list for exome-level data
          if ((class(sk) != "try-error") && (class(si) != "try-error")){
            class(sk) <- "skatCohort"
            message(log.tag, "Loaded Cohort ", cohort$name, ", SAMPLES=",
                    cohort$samples, ", GENES=", length(sk), ", SNPS=",
                    length(si$snpinfo$Name))
            skatcohort.exome.list <- lappend(skatcohort.exome.list, sk)
            snpinfo.exome.list <- lappend(snpinfo.exome.list, si)
            exome.cohorts.loaded <- c(exome.cohorts.loaded, cohort$name)
          }else{
            warning(log.tag,
                    "Error Loading Skatcohort Object Or Snpinfo File: COHORT=",
                    cohort$name, " CHR=", chrom, " TYPE=", cohort$type)
          }
        }else{
          warning(log.tag, "Missing Skatcohort Object Or Snpinfo File: COHORT=",
                  cohort$name, " CHR=", chrom, " TYPE=", cohort$type)
        }
      }
    }
    
    message(log.tag, "Exome Cohorts Loaded = ",
            paste(exome.cohorts.loaded, sep=",", collapse=","))
    
    #
    ## FOR EACH "GENE"-REGION ON CHROMOSOME
    #
    
    for (i in 1:length(regions.chr[,1])){
      cohorts.loaded <- exome.cohorts.loaded
      row <- regions.chr[i,]
      gene.start <- row$START
      gene.end <- row$END
      gene <- row$GENE_NAME
      log.tag <- paste("[",chrom, "_", gene,"] ", sep="")
      message(log.tag, "Processing Gene")
      ## Init data structs to collect info for gene-level skatMeta results
      snpinfo.list <- list()
      skatcohort.list <- list()

      ## Fetch gene-level data from previously fetched chromosome-level data
      
      ## Extract snp info for gene from chrom level objects
      for (item.id in 1:length(skatcohort.exome.list)){
        if (length(skatcohort.exome.list) == 0){
          message(log.tag, "There Appears To Be No Exome Data, Breaking From Loop")
          break;
        }
        snpinfo.item <- snpinfo.exome.list[[item.id]]
        skatcohort.item <- skatcohort.exome.list[[item.id]]
        t <- snpinfo.item$snpinfo[snpinfo.item$snpinfo$gene == gene, ]
        nsnps <- length(t[,1])
        message(log.tag, "Total Snps In Gene ", gene," For Cohort ",
                snpinfo.item$name, " = ", nsnps)
        if (nsnps > 0){
          snpinfo.list <- lappend(snpinfo.list, t)
        }else{
          warning(log.tag, "No Snps In Gene ", gene, " For Cohort ",
                  snpinfo.item$name, " = ", nsnps)
        }
        
        if ((length(skatcohort.item[[gene]]) == 5) &&
            (class(skatcohort.item[[gene]]$cov) == "matrix")){
          filtered.snps <- NULL
          sk <- list()
          sk[[gene]] <- skatcohort.item[[gene]]
          sk[[gene]]$cohortName <- snpinfo.item$name
          class(sk) <- "skatCohort"
          snpnames <- names(sk[[gene]]$maf)
          if (is.null(snpnames)){
            message(log.tag, "Exome: Failed to obtain SNP names from MAF attribute of skatCohort object, trying COV attribute.")
            snpnames <- colnames(sk[[gene]]$cov)
            if (is.null(snpnames)){
              warning(log.tag, "Exome: Failed to obtain SNP names from skatCohort object.")
            }
          }
          filtered.snps <- which(!(snpnames %in%  t$Name))
          colnames(sk[[gene]]$cov)[filtered.snps] <-
            paste(colnames(sk[[gene]]$cov)[filtered.snps], "x", sep="")
          rownames(sk[[gene]]$cov)[filtered.snps] <-
            paste(rownames(sk[[gene]]$cov)[filtered.snps], "x", sep="")
          if (!is.null(names(sk[[gene]]$maf))){
            names(sk[[gene]]$maf)[filtered.snps] <-
              paste(names(sk[[gene]]$maf)[filtered.snps], "x", sep="")
          }
          message(log.tag, "Filtered ", length(filtered.snps), " of ",
                  length(sk[[gene]]$maf),
                  " SNPs, remaining ", length(t$Name)," in SNPInfo")
          skatcohort.list <- lappend(skatcohort.list, sk)
        }else{
          warning(log.tag, "Invalid Exome-Level Skatcohort Gene Object For Gene ", gene, " For Cohort ",
                  snpinfo.item$name, " = ", nsnps)
        }
      }
      
      ## DEPRECATED: Load entire chromsome into skatMeta .. *slow*!
      ## skatcohort.list = skatcohort.exome.list

      #
      ## LOAD GENE-LEVEL DATA, TAGGED AS "GENOME"
      #
      
      for (cohort in cohorts){
        ## If cohort data is provided by chromosome ("genome")
        if (cohort$type == "genome"){
           ## Check files exist
          if (file.exists(cohort$skatcohort.file(chrom, gene)) &&
              file.exists(cohort$snpinfo.file(chrom, gene))) {
            message(log.tag, "Loading Cohort ", cohort$name)
            ## Load skatCohort and SNPInfo data
            sk <- try(cohort$skatcohort.load(cohort$skatcohort.file(chrom, gene)))
            si <- try(cohort$snpinfo.load(cohort$snpinfo.file(chrom, gene)))
            
            ## Check loading worked, append to list for gene-level data
            if ((class(sk) != "try-error") && (class(si) != "try-error")){
              if ((length(sk[[gene]]) != 5) ||
                  (class(sk[[gene]]$cov) != "matrix")){
                warning(log.tag, "Invalid Genome-Level Skatcohort Gene Object For Gene ", gene, " For Cohort ",
                        cohort$name)
                next;
              }
              snpnames <- names(sk[[gene]]$maf)
              if (is.null(snpnames)){
                message(log.tag, "Genome: Failed to obtain SNP names from MAF attribute of skatCohort object, trying COV attribute.")
                snpnames <- colnames(sk[[gene]]$cov)
                if (is.null(snpnames)){
                  warning(log.tag, "Genome: Failed to obtain SNP names from skatCohort object.")
                }
              }
              #sk[[gene]]$maf[!(snpnames %in%  si$Name)] <- 1.1
              filtered.snps <- which(!(snpnames %in%  si$Name))
              colnames(sk[[gene]]$cov)[filtered.snps] <-
                paste(colnames(sk[[gene]]$cov)[filtered.snps], "x", sep="")
              rownames(sk[[gene]]$cov)[filtered.snps] <-
                paste(rownames(sk[[gene]]$cov)[filtered.snps], "x", sep="")
              if (!is.null(names(sk[[gene]]$maf))){
                names(sk[[gene]]$maf)[filtered.snps] <-
                  paste(names(sk[[gene]]$maf)[filtered.snps], "x", sep="")
              }
              message(log.tag, "Filtered ", length(filtered.snps), " Of ",
                      length(sk[[gene]]$maf),
                      " SNPs, remaining ", length(si$Name)," In SNPInfo")

              

              sk[[gene]]$cohortName <- cohort$name
              skatcohort.list <- lappend(skatcohort.list, sk)
              snpinfo.list <- lappend(snpinfo.list, si)
              cohorts.loaded <- c(cohorts.loaded, cohort$name)
              message(log.tag, "Loaded Cohort ", cohort$name, ", Samples=",
                      cohort$samples,
                      ", Genes=", length(sk), ", Snps=", length(si$Name))
            }else{
              warning(log.tag,
                      "Unable To Load Skatcohort Object Or Snpinfo File Cohort=",
                      cohort$name, "  Chr=", chrom, " Gene=",gene, " Type=",
                      cohort$type)
            }
          }else{
            warning(log.tag, "One Or More Files Missing For COHORT=",
                    cohort$name, " CHR=", chrom, " GENE=",gene, " TYPE=",
                    cohort$type, " RDS.PATH=", cohort$skatcohort.file(chrom, gene), " SNPINFO.PATH=", cohort$snpinfo.file(chrom, gene))
          }
        }
      }

      ## #TODO FIX CASE WHERE THIS FAILS FOR SINGLE COHORT, LOW PRIORITY AS IT IS USED FOR LOGGING ONLY
      #message(log.tag, "COHORT LIST = ", paste(sapply(skatcohort.list, function(x) { x[[1]]$cohortName }), sep=",", collapse=","))
      #message(log.tag, "COHORT SAMPLE SIZES = ", paste(sapply(skatcohort.list, function(x) { x[[1]]$n }), sep=",", collapse=","))
      #message(log.tag, "TOTAL SAMPLE SIZE = ", sum(sapply(skatcohort.list, function(x) { x[[1]]$n })))
      
      message(log.tag, "Total Cohorts Loaded = ", paste(cohorts.loaded, sep=",",
                                                        collapse=","))
      message(log.tag, "Total Count Of Cohorts Loaded = ", length(cohorts.loaded))
      
      ## Bind all SNPInfo into one list
      snpinfo.set <- snpinfo.list
      
      snpinfo.list <- do.call("rbind", snpinfo.list)
      if (do.remove.single.cohort.singletons){
        snpinfo.list <- remove.single.cohort.singletons(snpinfo.list,
                                                        skatcohort.list)
      }
      nsnps <- length(snpinfo.list[,1])
      message(log.tag, "Total Snps Across Cohorts = ", nsnps)
      if (nsnps == 0){
        message(log.tag, "!!! NO SNP FOUND !!!")
        next;
      }
      ## Unique SNPInfo
      snpinfo.list <- as.data.frame(unique(snpinfo.list))
      ## Sort by position
      snpinfo.list <- snpinfo.list[with(snpinfo.list, order(Name)), ]
      
      nsnps <- length(snpinfo.list[,1])
      message(log.tag, "Total Unique Snps Across Cohorts = ",
              nsnps)

      ## Windows from SNPInfo
      p.value <- 1
      if (do.positional.windows){
        snpinfo.window.list <- get.positional.windows(snpinfo.list, gene.start, gene.end,
                                                      positional.window.size, positional.window.step, gene, log.tag)
        message(log.tag, "POSITIONAL WINDOW COUNT = ", length(snpinfo.window.list))
      }
      
      if (!do.positional.windows){
        snpinfo.window.list <- get.snp.windows(snpinfo.list, window.size, window.step, nsnps, gene)
        message(log.tag, "SNP WINDOW COUNT = ", length(snpinfo.window.list))
      }
      
      # win.id <- 1
      # for (win.s in seq(1, nsnps, by=window.step)){
      
      ## WINDOW ANALYSIS
      for (sw in snpinfo.window.list){
        log.tag <- paste("[",chrom, "_", gene, "_", sw$window.id,"] ", sep="")
        message(log.tag, "Running window meta-analysis")
        skatmeta.params <- list()

        ## !!! RESTART HERE !!!
        snpinfo.all.windows <- rbind(snpinfo.all.windows, sw$snpinfo)
        ## append skatCohort objects from each cohort to parameters for skatMeta
        nsubjects <- 0
        ## return(skatcohort.list)
        for (skatcohort.item in skatcohort.list){
          nsubjects <- nsubjects + skatcohort.item[[1]]$n
          skatmeta.params <- lappend(skatmeta.params, skatcohort.item)
        }
        ## return(skatcohort.list)
        ## Add other parameters
        if (skat.method == "singlesnpMeta"){
          skatmeta.params <- c(skatmeta.params, cohortBetas=FALSE)
          
        }
        if (skat.method != "burdenMeta" && skat.method != "singlesnpMeta"){
          skatmeta.params <- c(skatmeta.params, method=skat.pvalue.method)
        }
        if (skat.method != "singlesnpMeta"){
          skatmeta.params <- c(skatmeta.params, mafRange=skat.mafRange)
        }
        if (skat.method == "skatMeta"){
          skatmeta.params <- c(skatmeta.params, wts=skat.meta.wts)
        }
        if (skat.method == "skatOMeta"){
          skatmeta.params <- c(skatmeta.params, skat.wts=skat.meta.wts)
        }
        
        skatmeta.params <- c(skatmeta.params, SNPInfo=list(sw$snpinfo))
        ## return(skatcohort.list)
        if (skat.method == "singlesnpMeta"){
          out <- try(suppressWarnings(do.call(skat.method, skatmeta.params)))
        } else {
          out <- try(do.call(skat.method, skatmeta.params))
        }
        if ((class(out) !=  "try-error") || ((class(out) == "data.frame") && (!is.na(out$p)))) {
          maf <- if (is.null(out$maf)) "NA" else out$maf
          ntotal <- if (is.null(out$ntotal)) "NA" else out$ntotal
          name <- if (is.null(out$Name)) "NA" else out$Name
          cmaf <- if (is.null(out$cmaf)) "NA" else out$cmaf
          nsnps <- if (is.null(out$nsnps)) "NA" else out$nsnps
          qmeta <- if (is.null(out$Qmeta)) "NA" else out$Qmeta
          rho <- if (is.null(out$rho)) "NA" else out$rho
          pmin <- if (is.null(out$pmin)) "NA" else out$pmin
          errflag <- if (is.null(out$pmin)) "NA" else out$errflag
          beta <- if (is.null(out$beta)) "NA" else out$beta
          se <- if (is.null(out$se)) "NA" else out$se
          cmafTotal <- if (is.null(out$cmafTotal)) "NA" else out$cmafTotal
          cmafUsed <- if (is.null(out$cmafUsed)) "NA" else out$cmafUsed
          nsnpsTotal <- if (is.null(out$nsnpsTotal)) "NA" else out$nsnpsTotal
          nsnpsUsed <- if (is.null(out$nsnpsUsed)) "NA" else out$nsnpsUsed
          res.row <- data.frame(region=sw$region, gene=out$gene,
                                parent_region=out$gene, p.value=out$p,
                                ntotal=ntotal, name=name, maf=maf,
                                qmeta=qmeta, rho=rho, pmin=pmin, errflag=errflag,
                                beta=beta, se=se, cmafTotal=cmafTotal,
                                cmafUsed=cmafUsed, nsnpsTotal=nsnpsTotal,
                                nsnpsUsed=nsnpsUsed,
                                cmaf=cmaf, nmiss=out$nmiss, nsnps=nsnps,
                                chrom=chrom, window.start=sw$window.start,
                                window.end=sw$window.end, window.id=sw$window.id,
                                window.span=(sw$window.end - sw$window.start),
                                window.nsnps=length(sw$snpinfo$Name),
                                dropped.item = NA, cohorts = NA, nsubjects=nsubjects,
                                stringsAsFactors=FALSE)
          if (skat.method != "singlesnpMeta") {
            message(log.tag, paste(res.row, sep=" ",collapse=" "))
          }
          meta.res <- rbind(meta.res, res.row)
          p.value <- res.row$p.value
        }else{
          warning(log.tag, "SKATMETA ERROR OR INVALID P-VALUE")
          p.value <- 1
        }
        
        if (do.dropone && (skat.method != "singlesnpMeta") && (p.value <= dropone.max.pvalue)){
          message(log.tag, "SNP Drop-one On Window = ", sw$region)
          drop.id <- 1
          dropped.item <- ""
          ## For each SNP in window
          for (row.id in 1:length(sw$snpinfo$Name)){
            drop.id <- row.id
            snpinfo.window.drop <- sw$snpinfo[-row.id,]
            drop.snp <- sw$snpinfo[row.id,]
            snp.cohort <- ""
            cohort.str <- ""
            for (cohort.id in 1:length(snpinfo.set)){
              si <- snpinfo.set[[cohort.id]]
              if (drop.snp$Name %in% si$Name){
                snp.cohort <- paste(snp.cohort, "+", sep="")
                cohort.str <- paste(cohort.str,
                                    skatcohort.list[[cohort.id]][[1]]$cohortName,
                                    sep=",")
              }else{
                snp.cohort <- paste(snp.cohort, ".", sep="")
              }
            }
            snp.cohort <- paste(snp.cohort, cohort.str, sep=",")
            region.drop <- paste(gene, ".", sw$window.start, "-", sw$window.end, ".ds",
                                 drop.id, sep="")
            skatmeta.params[["SNPInfo"]] <- snpinfo.window.drop
            out <- try(do.call(skat.method, skatmeta.params))
            if ((class(out) !=  "try-error") || ((class(out) == "data.frame") && (!is.na(out$p)))) {
              ## if (class(out) !=  "try-error") {
              maf <- if (is.null(out$maf)) "NA" else out$maf
              ntotal <- if (is.null(out$ntotal)) "NA" else out$ntotal
              name <- if (is.null(out$Name)) "NA" else out$Name
              cmaf <- if (is.null(out$cmaf)) "NA" else out$cmaf
              nsnps <- if (is.null(out$nsnps)) "NA" else out$nsnps
              qmeta <- if (is.null(out$Qmeta)) "NA" else out$Qmeta
              rho <- if (is.null(out$rho)) "NA" else out$rho
              pmin <- if (is.null(out$pmin)) "NA" else out$pmin
              errflag <- if (is.null(out$pmin)) "NA" else out$errflag
              beta <- if (is.null(out$beta)) "NA" else out$beta
              se <- if (is.null(out$se)) "NA" else out$se
              cmafTotal <- if (is.null(out$cmafTotal)) "NA" else out$cmafTotal
              cmafUsed <- if (is.null(out$cmafUsed)) "NA" else out$cmafUsed
              nsnpsTotal <- if (is.null(out$nsnpsTotal)) "NA" else out$nsnpsTotal
              nsnpsUsed <- if (is.null(out$nsnpsUsed)) "NA" else out$nsnpsUsed
              res.row <- data.frame(region=region.drop, gene=out$gene,
                                    parent_region=sw$region, p.value=out$p,
                                    ntotal=ntotal, name=name, maf=maf,
                                    qmeta=qmeta, rho=rho, pmin=pmin,
                                    beta=beta, se=se, cmafTotal=cmafTotal,
                                    cmafUsed=cmafUsed,
                                    nsnpsTotal=nsnpsTotal, nsnpsUsed=nsnpsUsed,
                                    errflag=errflag, cmaf=cmaf,
                                    nmiss=out$nmiss, nsnps=nsnps,
                                    chrom=chrom, window.start=sw$window.start,
                                    window.end=sw$window.end, window.id=sw$window.id,
                                    window.span=(sw$window.end - sw$window.start),
                                    window.nsnps=length(sw$snpinfo$Name), 
                                    dropped.item = drop.snp$Name,
                                    cohorts = snp.cohort, nsubjects=nsubjects,
                                    stringsAsFactors=FALSE)
              # message(log.tag, paste(res.row, sep=" ",collapse=" "))
              meta.res <- rbind(meta.res, res.row)
            }else{
              warning(log.tag, "SKATMETA ERROR WHILE DROP SNP")
            }
            
          }  ## FOR

          #
          ## DROP ONE COHORT
          #
          
          drop.id <- 1
          message(log.tag, "Cohort Drop-one On Window = ", sw$region)
          for (row.id in 1:length(skatcohort.list)){
            skatmeta.params <- list()
            cohort.sym <- ""
            nsubjects <- 0
            dropped.item <- "NA"
            for (i in 1:length(skatcohort.list)){
              skatcohort.item <- skatcohort.list[[i]]
              if (i != row.id){
                skatmeta.params <- lappend(skatmeta.params, skatcohort.item)
                cohort.sym <- paste(cohort.sym, "+", sep="")
              }else{
                dropped.item <- skatcohort.item[[1]]$cohortName
                nsubjects <- skatcohort.item[[1]]$n
                cohort.sym <- paste(cohort.sym, ".", sep="")
              }
            }
            ## Add other parameters
            if (skat.method == "singlesnpMeta"){
              skatmeta.params <- c(skatmeta.params, cohortBetas=FALSE)
              
            }
            if (skat.method != "burdenMeta" && skat.method != "singlesnpMeta"){
              skatmeta.params <- c(skatmeta.params, method=skat.pvalue.method)
            }
            if (skat.method != "singlesnpMeta"){
              skatmeta.params <- c(skatmeta.params, mafRange=skat.mafRange)
            }
            if (skat.method == "skatMeta"){
              skatmeta.params <- c(skatmeta.params, wts=skat.meta.wts)
            }
            if (skat.method == "skatOMeta"){
              skatmeta.params <- c(skatmeta.params, skat.wts=skat.meta.wts)
            }

            
            skatmeta.params <- c(skatmeta.params, SNPInfo=list(sw$snpinfo))
            out <- try(do.call(skat.method, skatmeta.params))
            region.drop <- paste(gene, ".", sw$window.start, "-", sw$window.end, ".dc", drop.id,
                                 sep="")
            # browser()
            ##if (class(out) !=  "try-error") {
            if ((class(out) !=  "try-error") || ((class(out) == "data.frame") && (!is.na(out$p)))) {
              #singlesnp
              maf <- if (is.null(out$maf)) "NA" else out$maf
              ntotal <- if (is.null(out$ntotal)) "NA" else out$ntotal
              name <- if (is.null(out$Name)) "NA" else out$Name
              #skat
              cmaf <- if (is.null(out$cmaf)) "NA" else out$cmaf
              nsnps <- if (is.null(out$nsnps)) "NA" else out$nsnps
              qmeta <- if (is.null(out$Qmeta)) "NA" else out$Qmeta
              #skatO
              rho <- if (is.null(out$rho)) "NA" else out$rho
              pmin <- if (is.null(out$pmin)) "NA" else out$pmin
              errflag <- if (is.null(out$pmin)) "NA" else out$errflag
              #burden
              beta <- if (is.null(out$beta)) "NA" else out$beta
              se <- if (is.null(out$se)) "NA" else out$se
              cmafTotal <- if (is.null(out$cmafTotal)) "NA" else out$cmafTotal
              cmafUsed <- if (is.null(out$cmafUsed)) "NA" else out$cmafUsed
              nsnpsTotal <- if (is.null(out$nsnpsTotal)) "NA" else out$nsnpsTotal
              nsnpsUsed <- if (is.null(out$nsnpsUsed)) "NA" else out$nsnpsUsed

              res.row <- data.frame(region=region.drop, gene=out$gene,
                                    parent_region=sw$region, p.value=out$p,
                                    ntotal=ntotal, name=name, maf=maf,
                                    qmeta=qmeta, rho=rho, pmin=pmin,
                                    beta=beta, se=se, cmafTotal=cmafTotal,
                                    cmafUsed=cmafUsed,
                                    nsnpsTotal=nsnpsTotal, nsnpsUsed=nsnpsUsed,
                                    errflag=errflag, cmaf=cmaf,
                                    nmiss=out$nmiss, nsnps=nsnps,
                                    chrom=chrom, window.start=sw$window.start,
                                    window.end=sw$window.end, window.id=sw$window.id,
                                    window.span=(sw$window.end - sw$window.start),
                                    window.nsnps=length(sw$snpinfo$Name),
                                    dropped.item = dropped.item,
                                    cohorts = cohort.sym, nsubjects=nsubjects,
                                    stringsAsFactors=FALSE)
              
              # message(log.tag, paste(res.row, sep=" ",collapse=" "))
              meta.res <- rbind(meta.res, res.row)
            }else{
              warning(log.tag, "SKATMETA ERROR WHILE DROP COHORT")
            }
            drop.id <- drop.id + 1
          }
        }
        # win.id <- win.id + 1
      }
    }
  }
  return(list(meta=meta.res, snpinfo=snpinfo.all.windows))
}
