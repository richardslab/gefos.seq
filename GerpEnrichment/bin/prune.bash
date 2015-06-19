jid=$1
for i in {1..22}; do 
    echo "/usr/local/plink.1.90alpha/plink --bfile ${i}.maf0.005 \
--indep-pairwise 100k 20k 0.2 \
--out ${i}.maf0.005" |
    qsub -V -cwd -q all.q -N pr.${jid}$i
done
