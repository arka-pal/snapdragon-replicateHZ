#!/bin/bash

##### Relate: Estimate population size for each chromosome jointly with branch lengths and mutation rates
##### author: Arka Pal
##### written: 05.09.2024
##### update: 05.09.2024

##### USAGE: sbatch -J <jobname> job-relate_EstimatePopSize-joint-chrom.sbatch.sh <options>
##### $1: chrom

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
#SBATCH --mem-per-cpu=10G
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


## Load modules and set PATH
## -----
PATH_TO_RELATE=~/_softwares/relate_v1.2.2
export PATH=$PATH:~/_softwares/relate_v1.2.2/bin/:~/_softwares/relate_v1.2.2/scripts/
module load plink
module load bcftools
cd ~/snap_hap_repHZ/relate

## Initiate variables
## -----
baseDIR=~/snap_hap_repHZ/relate
chrom=$1
popLabels=$baseDIR/AvePla.MY.n74.poplabels
mu="5.7e-9"
Ne=813388
cd $baseDIR/Chr${chrom/chr}

echo -e "\n\n"
echo Population File: $popLabels
echo mutation rate: $mu
echo Ne: $Ne
echo -e "\n\n"

## Estimate population size per chromsome jointly with branch lengths, etc
## -----
echo -e '\n\nJoint estimation of population size and branch length for each chromosome\n-----'
srun time ~/_softwares/relate_v1.2.2/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i rel_$chrom \
              -o rel_${chrom}_joint \
              -m $mu \
              --poplabels $popLabels \
              --pops_of_interest AveM,AveY,PlaM,PlaY \
              --noanc 0 \
              --threshold 0 \
              --years_per_gen 3 \
              --num_iter 5 \
              --seed 420 \
              --threads ${SLURM_CPUS_PER_TASK}