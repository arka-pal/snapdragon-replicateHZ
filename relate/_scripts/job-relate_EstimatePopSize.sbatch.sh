#!/bin/bash

##### Relate: Estimate population size 
##### author: Arka Pal
##### written: 27.08.2024
##### update: 29.08.2024

##### USAGE: sbatch -J <jobname> job-relate_EstimatePopSize.sbatch.sh <options>


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

vcf=~/snap_hap_repHZ/statphase/AvePla_FrYe/$chrom.AvePla.FrYe.sorted.statphased.vcf.gz
# vcf=/nfs/scistore18/bartogrp/apal/snap_hap_repHZ/relate_test/test2.vcf.gz
ancestral=~/snap_hap_repHZ/ancestral_alleles/allelePolarised_chrom/ancestral_$chrom.txt
polarisedVcf=~/snap_hap_repHZ/relate/$chrom/$chrom.polarised
popLabels=$baseDIR/AvePla.MY.n74.poplabels

mu="5.7e-9"
Ne=813388


echo -e "\n\n"
echo Population File: $popLabels
echo mutation rate: $mu
echo Ne: $Ne
echo -e "\n\n"

## Symbolic link to all files by chromosome
ln -s ../*/rel_chr?.anc.gz .
ln -s ../*/rel_chr?.mut.gz .


## Estimate population size
## -----
echo -e '\n\nEstimating population size\n-----'
srun time $PATH_TO_RELATE/bin/RelateCoalescentRate \
                --mode EstimatePopulationSize \
                -m $mu \
                --poplabels $popLabels \
                -i rel \
                -o rel \
                --years_per_gen 3 \
                --first_chr 1 \
                --last_chr 8 \
                --num_samples 5 \
                --seed 420
