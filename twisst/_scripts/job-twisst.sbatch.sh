#!/bin/bash

##### Twisst analysis
##### author: Arka Pal
##### written: 22.07.2024
##### update: 05.09.2024

##### USAGE: sbatch -J <jobname> job-twisst.sbatch.sh <options>
# trees=$1
# outFile=$2
# popFile=$3

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
#SBATCH --mem-per-cpu=16G
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
export PATH=$PATH:$HOME/twisst/
module load python

baseDIR=~/snap_hap_repHZ/twisst/
trees=$1
outFile=$2
popFile=$3


## Run TWISST
## -----
srun time python ~/twisst/twisst.py -t $trees \
	-w $outFile \
	-g AveM -g AveY -g PlaM -g PlaY \
	--method fixed \
	--groupsFile $popFile