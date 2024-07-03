#!/bin/bash

##### SLURM script for bcftools mpileup & call
##### author: Arka Pal
##### written: 19.06.2024
##### update: 19.06.2024

##### Usage: sbatch -J <job-name> SnpCalling.sbatch.sh <options>
##### $1: baseDIR; ~/snap_hap_repHZ/SnpCalls
##### $2: bamlist; ~/snap_hap_repHZ/SnpCalls/
##### $3: outVCF


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
#SBATCH --cpus-per-task=5
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
baseDIR=$1
bamlist=$2
outVCF=$3
threads=${SLURM_CPUS_PER_TASK}
# baseDIR=~/snap_hap_repHZ/SnpCalls
refGenome=~/snap_hap/ref_genome/v3.5/Amajus_v3.5.fa

echo -e ref genome: $refGenome
echo -e bamFile: $bamlist
echo -e outVCF: $outVCF


## Run bcftools mpileup and call variants, and index
## -----
echo -e Calling SNPs
srun time   bcftools mpileup --threads $threads -b $bamlist \
                --annotate AD,ADF,ADR,DP,QS,SP -f $refGenome -Ou -d 500  | \
            bcftools call --threads $threads -Oz -o $outVCF \
                --annotate GQ,GP -m --write-index
## -----