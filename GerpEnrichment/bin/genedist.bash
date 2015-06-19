# Convert GWAS w/ GERP to BED format
tail -n +2 -q ../../merge/fn2stu.out.gerp | awk '{ print "chr"$2"\t"$1-1"\t"$1"\t"$3"\t"$20 }' > fn2stu.gerp.bed
# Use bedtools, get closest gene to GWAS SNP
/usr/local/bedtools-2.17.0/bin/bedtools closest -a fn2stu.gerp.bed -b ../knownGene.bed -d -t first > fn2stu.gerp.gene.bed

# Trim output
cut -f 1,2,3,5,10 fn2stu.gerp.gene.bed > fn2stu.gerp.gene.bed.trim

perl -p -e "s/\t/ /g;" fn2stu.gerp.gene.bed.trim > fn2stu.gerp.gene.bed.trim.space

paste -d ' ' fn2stu.gerp.gene.bed.trim.space <(tail -n +2 ../../merge/fn2stu.out.gerp) > fn2stu.gerp.gene.bed.trim.space.meta

cut -d ' ' -f 1,2,3,4,5,8,11,17 fn2stu.gerp.gene.bed.trim.space.meta > fn2stu.gerp.gene.bed.trim.space.meta.trim


