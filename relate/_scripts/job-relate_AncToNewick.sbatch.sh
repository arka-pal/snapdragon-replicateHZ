#!/bin/bash

##### Convert Relate trees to newick
##### author: Arka Pal
##### written: 05.09.2024
##### update: 05.09.2024

##### USAGE: sbatch -J <jobname> job-relate-AncToNewick.sbatch.sh <options>
# chrom=$1
# start=$2
# end=$3
# ancFile=$4
# mutFile=$5
# outPrefix=$6


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


## Initiate variables
## -----
baseDIR=~/snap_hap_repHZ/relate
chrom=$1
start=$2
end=$3
ancFile=$4
mutFile=$5
outPrefix=$6

echo -e "\n\n"
echo chrom: $chrom
echo "start:" $start
echo "end:" $end
echo .anc file: $ancFile
echo .mut file: $mutFile
echo out prefix: $outPrefix
echo -e "\n\n"


## Convert to newick trees
## -----
srun $PATH_TO_RELATE/bin/RelateExtract\
                 --mode AncToNewick \
                 --anc $ancFile \
                 --mut $mutFile \
                 --first_bp $start \
                 --last_bp $end \
                 -o $outPrefix

## Clean up
## -----
bgzip $outPrefix.newick