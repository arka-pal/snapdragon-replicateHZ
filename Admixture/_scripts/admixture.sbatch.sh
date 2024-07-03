#!/bin/bash

##### Admixture Analysis
##### author: Arka Pal
##### written: 27.06.2024
##### update: 27.06.2024

##### USAGE: sbatch -J <jobname> job-admixture.sbatch.sh <options>
##### $1 - inFile
##### $2 - K


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------
##### NODE ISSUES to check
##### SBATCH --exclude=zeta[243-263],bigterra152,beta[231-235]
##### SBATCH --gres=gpu:1
##### -----



## Load module and set path
## -----
# module load plink
export PATH=$PATH:~/_softwares/admixture_linux-1.3.0


## Read variables
## -----
input=$1
K=$2


## Run admixture
## -----
srun admixture $input $K
