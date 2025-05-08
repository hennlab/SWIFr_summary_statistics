# SWIFr summary statistics computation
Pipeline to compute genome-wide summary statistics to detect balancing selection

This pipeline runs on a VCF with 1 chromosome and computes summary stats on window sizes of 100kb.
The script runs as a job array, splitting the VCF into regions of 8Mb and computing the summary statistics in parallel to optimoze computation time 
Output is a CSV file with summary stats per variant position (user-defined minimum variant frequency)

Summary statistics based on BaSe ([Isildak, et al 2021](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13379))

## Steps

### 1.  Define variables in bash
```
chr= #chromosome. Format: integer \
path="" #path to VCF infput. Format: string \
vcf="" #name of VCF. Format: string \
out="" #output name. Format: string \
freq= #minimum frequency for a variant to be considered a candidate (aka to calculate summ. stats on). Format: float
```

### 2. Create bed file with ranges of 8Mb to run in a job array (start and end regions separated by whitespace)
```
module load bcftools; bcftools query -f '%POS' $path/$vcf -o chr${chr}_position.txt
./create_bed.sh chr${chr}_position.txt chr${chr}_positions.bed
rm chr${chr}_position.txt
```

### 4. Run runstats script to launch the jobs
```
mkdir logs/
num=$(wc -l chr${chr}_positions.bed | cut -d' ' -f1) #number of arrays depending on number of regions to split vcf
sbatch --array=${num}-${num}%20 1.runstats_array.sh $chr $path $vcf $freq $out
```

### 5. Track progress of job array: 
```
#1. count number of variants done
grep "Processing" logs/OUT_runstats_100kb_* | wc -l 
#2. count number of successful slurm arrays done
tail -n 1 logs/OUT_runstats_100kb_* | grep "Successful" | wc -l 
#3. count number of non-successful (or running) slurm arrays
tail -n 1 logs/OUT_runstats_100kb_* | grep -v 'Successful\|=\|^$' | wc -l
```

### 6. Merge output from each slurm array into one CSV file
```
head -n1 ${out}.chr${chr}.1.mod.csv > ${out}.chr${chr}.all.csv
tail -n+2 ${out}.chr${chr}.*.mod.csv | grep -v "==>" >> ${out}.chr${chr}.all.csv
```
