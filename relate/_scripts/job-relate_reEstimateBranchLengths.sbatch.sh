#!/bin/bash

##### Relate: ReEstimate branch lengths
##### author: Arka Pal
##### written: 05.09.2024
##### update: 05.09.2024

##### USAGE: sbatch -J <jobname> job-relate_reEstimateBranchLengths.sbatch.sh <options>
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
inPrefix=$baseDIR/reEstimateBranchLengths/rel_$chrom
outPrefix=$baseDIR/reEstimateBranchLengths/rel_${chrom}_updated
mrate=$baseDIR/estimateMutRate/rel_avg.rate
coal=$baseDIR/estimatePopSize_joint/rel_joint.coal
mu="5.7e-9"

echo -e "\n\n"
echo inPrefix: $inPrefix
echo outPrefix: $outPrefix
echo mrate: $mrate
echo coal: $coal
echo mutation rate: $mu
echo -e "\n\n"

# ## Symbolic link to all files by chromosome
# ln -s ../*/rel_chr?.anc.gz .
# ln -s ../*/rel_chr?.mut.gz .


## Reestimate branch lengths
## -----
echo -e '\n\nreEstimating branch length\n-----'
srun time $PATH_TO_RELATE/bin/RelateCoalescentRate \
                --mode ReEstimateBranchLengths \
                -i $inPrefix \
                -o $outPrefix \
                --mrate $mrate \
                --coal $coal \
                -m $mu \
                --seed 420

## Clean up
## -----
bgzip rel_${chrom}_updated*