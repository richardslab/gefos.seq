#!/bin/bash

mkdir -p out/
mkdir -p log/
echo 'Rscript run_meta.devel.R ~/common/gefos.seq/SKATMETA/static/batched_regions_autosomes/${SGE_TASK_ID}.txt cohorts.R > out/${SGE_TASK_ID}.txt' | qsub -t 1-248 -N $1 -cwd -V -q all.q -o log -e log