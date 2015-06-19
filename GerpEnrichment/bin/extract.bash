jid=$1
for i in {1..22}; do 
    echo "/usr/local/plink.1.90alpha/plink --bfile ${i}.bi \
    --extract snps.txt \
    --out ${i}.maf0.005 \
    --make-bed" |
    qsub -V -cwd -q all.q -N j${jid}.$i
done

