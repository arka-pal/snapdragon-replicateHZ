#!/bin/bash

##### Twisst analysis in windows
##### author: Arka Pal
##### written: 23.07.2024
##### update: 23.07.2024

##### USAGE: sbatch -J <jobname> job-twisst_windows.sbatch.sh <options>


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%a_%A.out
#SBATCH --error=%x_%a_%A.out
#SBATCH --open-mode=append
#SBATCH --partition=gpu
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
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


## Read input and set variables
## -----
baseDIR=~/snap_hap_repHZ/twisst
chrom=$1
treesALL=$2
popFile=$3

windowID=${SLURM_ARRAY_TASK_ID}
windowSize=10000
windowStart=$((($windowID-1)*$windowSize+1))
windowEnd=$(($windowStart+$windowSize-1))

treesWindow=$baseDIR/trees_tskit/trees_in_windows/$chrom.$windowID.$windowStart.$windowEnd.AvePla.all.treeList.txt
outputWindow=$baseDIR/trees_tskit/twisst_in_windows/$chrom.$windowID.$windowStart.$windowEnd.weights.csv.gz



## Create trees in windows
echo -e $windowID'\t'$windowStart'\t'$windowEnd'\n\n'
sed -n "${windowStart},${windowEnd}p" $treesALL > $treesWindow


## Run TWISST
## -----
python ~/twisst/twisst.py -t $treesWindow \
	-w $outputWindow \
	-g Ave_FR -g Ave_Y -g Pla_FR -g Pla_Y \
	--method fixed \
	--groupsFile $popFile

rm $treesWindow