jid=$1
for i in {1..22}; do 
    echo "plink --bfile ${i}.maf0.005 \
--noweb --allow-no-sex --r2 --ld-window 1000 --ld-window-kb 1000 --nonfounders \
--ld-window-r2 0.2 --ld-snp-list ${i}.maf0.005.prune.in --out ${i}.maf0.005" \
	| qsub -V -cwd -q all.q -N r2.${jid}.$i
done
