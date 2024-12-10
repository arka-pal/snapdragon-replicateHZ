#!/bin/bash

##### Relate: Run singer on windows
##### author: Arka Pal
##### written: 03.09.2024
##### update: 05.09.2024

##### USAGE: sbatch -J <jobname> job-singer_multiWindows.sbatch.sh <options>


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%a_%A.out
#SBATCH --error=%x_%a_%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
#SBATCH --time=12:00:00
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
# module load singer
module load vcftools
module load bcftools

## Initiate variables
## -----
baseDIR=~/snap_hap_repHZ/singer
chrom=$1
# chrom=Chr2
winFile=$2
# winFile=$baseDIR/input/win100k_SHM_Ne.tsv

windowID=${SLURM_ARRAY_TASK_ID}
start=$(sed -n "$((windowID+2))p" $winFile | cut -f2)
end=$(sed -n "$((windowID+2))p" $winFile | cut -f3)

inVcf=$baseDIR/input/Flavia.sg
outPrefix=$baseDIR/runs/$chrom.win$windowID.$start.$end/out
if [ ! -d $baseDIR/runs/$chrom.win$windowID.$start.$end ]; then mkdir $baseDIR/runs/$chrom.win$windowID.$start.$end; fi

mu=5.7e-9
ratio=1
Ne=$(sed -n "$((windowID+2))p" $winFile | cut -f6) #vcftools
# Ne=$(sed -n "$((windowID+2))p" $winFile | cut -f8) #SHM

mcmc_iters=100
thin=20
polar=0.5

echo "VCF: $inVcf.vcf"
echo "chrom: $chrom"
echo "window: $windowID"
echo "start: $start"
echo "end: $end"
echo "outPrefix: $outPrefix"
echo "mu: $mu"
echo "Ne: $Ne"
echo "mcmc iters: $mcmc_iters"
echo "mcmc thin: $thin"
echo "polar: $polar"
echo -e "-----\n\n\n"


## Run singer 
#NB: ~150mins on head node; 500KB windows
srun time ~/_softwares/SINGER/releases/singer-0.1.7-beta-linux-x86_64/singer_master -vcf $inVcf \
                    -output $outPrefix \
                    -start $start \
                    -end $end \
                    -m $mu \
                    -ratio $ratio \
                    -Ne $Ne \
                    -n $mcmc_iters \
                    -thin $thin \
                    -polar $polar \
                    -seed 420
