#!/bin/bash

##### Infer trees with tsinfer
##### author: Arka Pal
##### written: 22.07.2024
##### update: 22.07.2024

##### USAGE: sbatch -J <jobname> job-ts.sbatch.sh <options>
##### $1
##### $2
##### $3
##### $4
##### $5
##### $6


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
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
module load python
python ~/snap_hap_repHZ/tsinfer/_scripts/run_ts.py $1 $2 $3 $4 $5 $6