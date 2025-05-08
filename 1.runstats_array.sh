#!/bin/bash
#SBATCH -n 1
#SBATCH --job-name=stats
#SBATCH -t 4-00:00:00 # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH -o logs/OUT_runstats_100kb_%a.out
#SBATCH --mem=12G

source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate tskit
chr=$1
path=$2
vcf=$3
freq1=$4
out_1=$5

## Find regions to split vcf based on ${SLURM_ARRAY_TASK_ID}
var=$(sed -n -e ${SLURM_ARRAY_TASK_ID}p chr${chr}_positions.bed)
start=$(echo $var | sed 's/ .*$//g' )
end=$(echo $var | sed 's/^.* //g' )
last=$(wc -l chr${chr}_positions.bed | cut -d' ' -f1)

samplesize=$(grep "CHROM" ${path}/${vcf} | cut -f10- | tr '\t' '\n' | wc -l)

# Outpuf file name
out="${out_1}.chr${chr}.${SLURM_ARRAY_TASK_ID}.csv" 

# Initialize return code
rc=0

## Run stats script
python -u 2.stats_from_vcf_CAAPA_freq_array.py $path $vcf $start $end $samplesize $freq1 $out ${SLURM_ARRAY_TASK_ID} $last

# Capture the exit status
rc=$?

# Check the exit status and print the appropriate message
if [[ $rc -eq 0 ]]; then
    echo "Successful run"
else
    echo "Run failed with exit code $rc"
fi
