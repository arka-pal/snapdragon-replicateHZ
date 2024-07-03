#!/bin/bash

##### SLURM script to postprocess STITCH output 
##### author: Arka Pal
##### written: 25.06.2024
##### update: 25.06.2024

##### USAGE: sbatch -J <jobname> PostStitch.sbatch.sh <options>
##### $1 - baseDIR
##### $2 - chrom


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
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



## Set path
## -----
export PATH=~/.local/bin:$PATH


## Read and define variables
## -----
baseDIR=$1
chrom=$2
stitchVcfList=$baseDIR/$chrom/stitchVCFs.highConf.list
# ncores=${SLURM_CPUS_PER_TASK}
ncores=4
cd $baseDIR/$chrom


## Postprocessing
## -----
# Make list of STITCH VCFs
realpath $baseDIR/$chrom/stitch_chromSegments/*/stitch.*.vcf.gz  > $stitchVcfList

# Index all STITCHVCFs
for vcf in `cat $stitchVcfList`; do echo $vcf; bcftools tabix -f $vcf; done

# Concatenate all STITCH VCFs; index; add AC+AN tags and extract STITCH-parameters
highConfVcf=$baseDIR/$chrom/$chrom.AvePla.stitch.vcf.gz
bcftools concat --allow-overlaps --threads $ncores -Oz -o $highConfVcf -f $stitchVcfList
bash ~/snap_hap/_scripts/bash/fill-AC-AN_vcfgzFormat.sh $highConfVcf
mv ${highConfVcf/.vcf.gz/.tagged.vcf.gz} $highConfVcf
bcftools tabix -f $highConfVcf
bcftools query -f "%CHROM\t%POS\t%EAF\t%INFO_SCORE\t%HWE\t%ERC\t%EAC\t%PAF\t%AN\t%AC" $highConfVcf > ${highConfVcf/.vcf.gz/.params}

# Extract all genotypes
bash ~/snap_hap/_scripts/bash/stitch/add_PG_PL.sh $highConfVcf

# Add AC+AN tags and extract STITCH-parameters
FullVcf=$baseDIR/$chrom/$chrom.AvePla.stitch.PL.vcf.gz
bash ~/snap_hap/_scripts/bash/fill-AC-AN_vcfgzFormat.sh $FullVcf
mv ${FullVcf/.vcf.gz/.tagged.vcf.gz} $FullVcf
bcftools tabix -f $FullVcf
bcftools query -f "%CHROM\t%POS\t%EAF\t%INFO_SCORE\t%HWE\t%ERC\t%EAC\t%PAF\t%AN\t%AC" $FullVcf > ${FullVcf/.vcf.gz/.params}

# Keep only variant sites
FinalVcf=$baseDIR/$chrom/$chrom.AvePla.stitch.SnpOnly.final.vcf.gz
bcftools view $FullVcf -e 'AC==0 | AC==AN' -Oz -o $FinalVcf --threads $ncores
bcftools tabix -f $FinalVcf
bcftools query -f "%CHROM\t%POS\t%EAF\t%INFO_SCORE\t%HWE\t%ERC\t%EAC\t%PAF\t%AN\t%AC" $FinalVcf > ${FinalVcf/.vcf.gz/.params}

rm *header
rm *tagged.vcf.gz.tbi