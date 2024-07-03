#!/bin/bash

##### Genome scan with genomics_general
##### author: Arka Pal
##### written: 21.06.2024
##### update: 21.06.2024

##### Usage: sbatch -J <job-name> WindowStats_gen.sbatch.sh <options>
##### $1: geno
##### $2: output



## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
### #SBATCH --exclude=zeta[243-263],beta[231-235]

#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G

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
module load bcftools
module load python
export PATH=$PATH:$HOME/genomics_general:$HOME/genomics_general/VCF_processing


## Read input & define variables
## -----
# geno=./test.geno.gz
# region=Chr6
# # output=PlaFR-PlaY_Chr6_coords_win3kb_step3kb.csv.gz
# output=test.csv.gz
# window=3000
# pool1=Pla_FR
# pool2=Pla_Y
# popFile=./AvePla_Pools.txt


## Estimate pi, dxy and Fst
## -----
srun python ~/genomics_general/popgenWindows.py -g ~/snap_hap_repHZ/genome_scans/geno/AvePla.Chr6.biSNPs.filtered.test.geno.gz \
	-o ~/snap_hap_repHZ/genome_scans/PlaFR-PlaY_Chr6_coords_win10kb_step10kb.csv.gz \
	-f diplo \
	--windType coordinate \
    --windSize 10000 \
    --stepSize 10000 \
	-p Pla_FR -p Pla_Y \
	--popsFile ~/snap_hap_repHZ/genome_scans/AvePla-FRYe_pools.txt \
	--writeFailedWindow \
    --addWindowID \
    --threads 5
echo -e "\nDone"
## -----