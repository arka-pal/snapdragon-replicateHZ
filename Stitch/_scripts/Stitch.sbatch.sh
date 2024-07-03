#!/bin/bash

##### SLURM script to run STITCH imputation 
##### author: Arka Pal
##### written: 23.06.2024
##### update: 23.06.2024

##### USAGE: sbatch --array=1-n -J <job-name> Stitch.sbatch.sh <options>
##### $1 - chrom
##### $2 - buffer
##### $3 - K #75
##### $4 - downsampleToCov #20
##### $5 - use_bx_tag (TRUE/FALSE) #TRUE
##### $6 - ngen #100
##### $7 - niter #40
##### $8 - expRate #0.5
##### $9 - plot (TRUE/FALSE) #TRUE


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%a-%A.out
#SBATCH --error=%x_%a-%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=16G

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
module load stitch bcftools


## Read and define variables
## -----
baseDIR=~/snap_hap_repHZ
bamlist=$baseDIR/SnpCalls/AvePla.bams.list

# STITCH regions
chrom=$1
chromSegments=~/snap_hap/ref_genome/chromSegments/${chrom}_segments.txt
window=${SLURM_ARRAY_TASK_ID}
chromStart=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $chromSegments | cut -f1)
chromEnd=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $chromSegments | cut -f2)
buffer=$2

# STITCH variables
posfile=$baseDIR/Stitch/$chrom/${chrom}.AvePla.biSNPs.filtered.missLT80.pos
K=$3
downsampleToCov=$4
use_bx_tag=$5
ngen=$6
niter=$7
expRate=$8
plot=$9
outputDIR=$baseDIR/Stitch/$chrom/stitch_chromSegments/${chrom}_${chromStart}-${chromEnd}_buffer$buffer_K${K}_cov${downsampleToCov}_bxTRUE_niter${niter}_ngen${ngen}_r${expRate}_plotFALSE


## Run STITCH
## -----
echo bamlist: $bamlist
echo posfile: $posfile
echo OUTPUT: $outputDIR
echo Region: $chrom $chromStart $chromEnd $buffer

srun STITCH.R --chr=$chrom \
    --K=$K \
    --downsampleToCov=$downsampleToCov \
    --use_bx_tag=${use_bx_tag} \
    --niterations=$niter \
    --expRate=$expRate  \
    --nGen=$ngen \
    --plotAfterImputation=$plot \
    --plotHapSumDuringIterations=$plot \
    --plot_shuffle_haplotype_attempts=$plot \
    --regionStart=$chromStart \
    --regionEnd=$chromEnd \
    --buffer=$buffer \
    --bamlist=$bamlist \
    --posfile=$posfile \
    --outputdir=$outputDIR \
    --nCores=${SLURM_CPUS_PER_TASK}

