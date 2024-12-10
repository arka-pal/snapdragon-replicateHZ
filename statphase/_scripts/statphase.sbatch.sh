#!/bin/bash

##### SLURM script: Statistical Phasing with shapeit5
##### author: Arka Pal
##### written: 01.07.2024
##### update: 01.07.2024

##### USAGE: sbatch -J <job-name> statphase.sbatch.sh <options>
##### $1 - targetVcf
##### $2 - outVcf
##### $3 - chrom


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-263],bigterra152,beta[231-235]
#SBATCH --time=240:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END
#SBATCH --no-requeue
#SBATCH --export=NONE
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------
##### NB: NODE ISSUES; check if to be included -----
##### SBATCH --exclude=zeta[243-263],bigterra152,beta[231-235]
##### SBATCH --gres=gpu:1
##### -----


## Load modules
## -----
module load shapeit5 bcftools vcftools


## Initiate variables
## -----
targetVcf=$1
outVcf=$2
chrom=$3


## Shapeit5 without ref panel
## -----
phase_common_static --input $targetVcf --output $outVcf --region $chrom --thread ${SLURM_CPUS_PER_TASK}