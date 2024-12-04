#!/bin/bash

##### Full script for Relate
##### author: Arka Pal
##### written: 26.08.2024
##### update: 26.08.2024

##### USAGE: sbatch -J <jobname> job-relate.sbatch.sh <options>
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
#SBATCH --mem-per-cpu=12G
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
export PATH=$PATH:~/_softwares/relate_v1.2.2/bin/:~/_softwares/relate_v1.2.2/scripts/
module load bcftools
module load plink

## Initiate variables
## -----
chrom=$1
baseDIR=~/snap_hap_repHZ/relate

vcf=~/snap_hap_repHZ/statphase/AvePla_FrYe/$chrom.AvePla.FrYe.sorted.statphased.vcf.gz
# vcf=/nfs/scistore18/bartog/apal/snap_hap_repHZ/relate_test/test2.vcf.gz
ancestral=~/snap_hap_repHZ/ancestral_alleles/allelePolarised_chrom/ancestral_$chrom.txt
polarisedVcf=~/snap_hap_repHZ/relate/$chrom/$chrom.polarised
popLabels=$baseDIR/AvePla.MY.n74.poplabels

mu="5.7e-9"
Ne=813388

if [ ! -d $baseDIR/$chrom ]; then mkdir -p $baseDIR/$chrom; fi
cd $baseDIR/$chrom

echo -e "\n\n"
echo VCF: $vcf
echo ancestral: $ancestral
echo Polarised VCF: $polarisedVcf
echo Population File: $popLabels
echo mutation rate: $mu
echo Ne: $Ne
echo -e "\n\n"

## Step 1. Create text file with ancestral allele information ##
## -----
echo -e '\n\nStep 1. Create text file with ancestral allele information\n-----'
tail +2 ~/snap_hap_repHZ/ancestral_alleles/allelePolarised_chrom/allelePolarised_${chrom}.csv | \
    cut -d, -f2,1,35 | tr ',' '\t' | awk '{print $2":"$1"\t"$3}' > $ancestral

## Step 2. Recode VCF to ensure ancestral allele is always coded by 0 ##
## -----
echo -e '\n\nStep 2. Recode VCF to ensure ancestral allele is always coded by 0\n-----'
plink2 --vcf $vcf --set-all-var-ids Chr@:\# --ref-allele 'force' $ancestral 2 1 --export vcf --out $polarisedVcf
bgzip -f $polarisedVcf.vcf

## Step 3. Convert VCF to haps/sample format ##
## -----
echo -e '\n\nStep 3. Convert VCF to haps/sample format\n-----'
RelateFileFormats --mode ConvertFromVcf --haps $chrom.haps --sample $chrom.sample -i $polarisedVcf

## Step 4. Generate SNP annotation ##
## -----
echo -e '\n\nStep 4. Generate SNP annotation\n-----'
RelateFileFormats --mode GenerateSNPAnnotations --haps $chrom.haps --sample $chrom.sample --poplabels $popLabels -o $chrom

## Step 5. Run relate
## NB: recombination map created manually.
#time Relate --mode All -m $mu -N $Ne --haps $chrom.haps --sample $chrom.sample --map $chrom.map --annot $chrom.annot --seed 420 -o $chrom
# Run Relate in parallel
## -----
echo -e '\n\nStep 5. Run relate\n-----'
srun time ~/_softwares/relate_v1.2.2/scripts/RelateParallel/RelateParallel.sh --mode All \
                            -m $mu -N $Ne \
                            --haps $chrom.haps \
                            --sample $chrom.sample \
                            --map $chrom.map \
                            --annot $chrom.annot \
                            --seed 420 \
                            -o rel_chr${chrom/Chr} \
                            --threads ${SLURM_CPUS_PER_TASK}

## Estimate population size
## -----
echo -e '\n\nEstimating population size\n-----'
srun time ~/_softwares/relate_v1.2.2/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i rel_chr${chrom/Chr} \
              -m $mu \
              --poplabels $popLabels \
              --pop_of_interest AveM,AveY,PlaM,PlaY \
              --noanc 0 \
              --seed 420 \
              --threshold 0 \
              -o rel_chr${chrom/Chr}.popsize \
              --threads ${SLURM_CPUS_PER_TASK}

# ## Convert Relate output to tskit format
# ## -----
# srun time ~/_softwares/relate_lib/bin/Convert --mode ConvertToTreeSequence \
#               --compress \
#               --anc $chrom.anc \
#               --mut $chrom.mut \
#               -o $chrom
