extract.bash -- extract SNPs with MAF > 0.5%
prune.bash -- Prune SNPs for LD
r2.bash -- Calculate R2 for pruned SNPs
r2clean.bash -- Fix tab/space issues in plink R2 output file
merge.awk -- Merge R2 and p-values
select.awk -- Select lowest P-value
filtermeta.awk -- Filter meta analysis SNPs for those  from select.awk (lowest independant p-value)
matchit.R -- Run matchit and FDR analysis
table.R -- generate table from FDR results
