#!/bin/bash

##### SLURM script for SNP filtering
##### author: Arka Pal
##### written: 20.06.2024
##### update: 23.06.2024

##### Usage: sbatch -J <job-name> SnpFiltering.sbatch.sh <options>
##### $1: inVCF
##### $2: outVCF


## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
### #SBATCH --exclude=zeta[243-263],beta[231-235]

#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------
##### NB: NODE ISSUES; check if to be included -----
##### SBATCH --exclude=zeta[243-263],bigterra152,beta[231-235]
##### SBATCH --gres=gpu:1
##### -----


## Load modules
## -----
module load bcftools


## Read input & define variables
## -----
inVCF=$1
outVCF=$2
threads=${SLURM_CPUS_PER_TASK}
# baseDIR=~/snap_hap_repHZ/SnpCalls

echo -e inVCF: $inVCF
echo -e outVCF: $outVCF


## Filter sites
## -----
echo -e "Filter Step 1: keep biSNPs only" 
srun time bcftools filter --threads $threads $inVCF --SnpGap 5 | \
     bcftools view --threads $threads - -Oz -o $outVCF -m2 -M2 -v snps -e "AC==0 || AC==AN" --write-index

echo -e "Filter Step 2: keep biSNPs based on QUAL, MQ, DP" 
srun time bcftools filter --threads $threads $inVCF -e "INFO/DP>130 | QUAL<20 | MQ<20" -Oz -o $outVCF --write-index

echo -e "Filter Step 3: Remove SNPs with missing fraction > 0.8" 
srun time bcftools filter --threads $threads $inVCF -e "F_MISSING>0.80" -Oz $outVCF --write-index
echo -e "\nDone"

echo -e "Make pos file" 
srun time bcftools view -H $outVCF | cut -f1,2,4,5 > ${outVCF/.vcf.gz/.pos}
echo -e '\nFinal no. of sites:' `wc -l ${outVCF/.vcf.gz.pos}`
echo -e "\nDone"

## -----