##### Draw trees from Relate

## Libraries and paths ----
setwd('~/snap_hap_repHZ/relate/')
library(foreach)
library(tidyverse)
library(data.table)
library(tictoc)
library(patchwork)
library(ape)
library(ggtree)
source('~/twisst/plot_twisst.R')
source('./_scripts/treeView.R')
## --


## Read smoothed weights ----
chr1 = fread('~/snap_hap_repHZ/twisst/trees_relate/weights/chr1.weights.smooth.csv.gz')
chr2 = fread('~/snap_hap_repHZ/twisst/trees_relate/weights/chr2.weights.smooth.csv.gz')
chr3 = fread('~/snap_hap_repHZ/twisst/trees_relate/weights/chr3.weights.smooth.csv.gz')
chr4 = fread('~/snap_hap_repHZ/twisst/trees_relate/weights/chr4.weights.smooth.csv.gz')
chr5 = fread('~/snap_hap_repHZ/twisst/trees_relate/weights/chr5.weights.smooth.csv.gz')
chr6 = fread('~/snap_hap_repHZ/twisst/trees_relate/weights/chr6.weights.smooth.csv.gz')
chr7 = fread('~/snap_hap_repHZ/twisst/trees_relate/weights/chr7.weights.smooth.csv.gz')
chr8 = fread('~/snap_hap_repHZ/twisst/trees_relate/weights/chr8.weights.smooth.csv.gz')
## --


## Read trees ----
# tr1 = fread('~/snap_hap_repHZ/relate/newickTrees/rel_chr1-1-71919034.newick.gz', sep='', header=FALSE)
# colnames(tr1) = 'newickTree'
tr2 = fread('~/snap_hap_repHZ/relate/newickTrees/rel_chr2-1-77118269.newick.gz', sep='', header=FALSE)
colnames(tr2) = 'newickTree'
# tr3 = fread('~/snap_hap_repHZ/relate/newickTrees/rel_chr3-1-65231163.newick.gz', sep='', header=FALSE)
# colnames(tr3) = 'newickTree'
tr4 = fread('~/snap_hap_repHZ/relate/newickTrees/rel_chr4-1-54887108.newick.gz', sep='', header=FALSE)
colnames(tr4) = 'newickTree'
tr5 = fread('~/snap_hap_repHZ/relate/newickTrees/rel_chr5-1-71106538.newick.gz', sep='', header=FALSE)
colnames(tr5) = 'newickTree'
tr6 = fread('~/snap_hap_repHZ/relate/newickTrees/rel_chr6-1-55699338.newick.gz', sep='', header=FALSE)
colnames(tr6) = 'newickTree'
# tr7 = fread('~/snap_hap_repHZ/relate/newickTrees/rel_chr7-1-55564713.newick.gz', sep='', header=FALSE)
# colnames(tr7) = 'newickTree'
# tr8 = fread('~/snap_hap_repHZ/relate/newickTrees/rel_chr8-1-57431585.newick.gz', sep='', header=FALSE)
# colnames(tr8) = 'newickTree'
## --


## combine trees and weights ----
# dat1 = cbind(chr1, tr1)
dat2 = cbind(chr2, tr2)
# dat3 = cbind(chr3, tr3)
dat4 = cbind(chr4, tr4)
dat5 = cbind(chr5, tr5)
dat6 = cbind(chr6, tr6)
# dat7 = cbind(chr7, tr7)
# dat8 = cbind(chr8, tr8)
## --


## Tip colors ----
tip.colors = c(rep('magenta',38), rep('yellow3',38), rep('pink3',36), rep('orange',36))
ids = 0:147
tipInfo = data.frame(IDs = ids, tipColors = tip.colors)
tipInfo
## --


## ROS1 trees ----
source('./_scripts/treeView.R')
geneName = 'ros' 
PATH_TO_RELATE = '~/_softwares/relate_v1.2.2'
filename_haps = './Chr6/Chr6.haps.gz'
filename_sample = './Chr6/Chr6.sample'
filename_poplabels = './AvePla.MY.n74.poplabels'
filename_anc = './Chr6/rel_chr6_joint.anc.gz'
filename_mut = './Chr6/rel_chr6_joint.mut.gz'
years_per_gen = 3
# snpPos = c(52884457, 52884489, 52884528, 52884553, 52884570, 52884624, 52884770)
treeIDs = 513855:513861
# snpPos = 52884457
rosList = list()
for (treeID in treeIDs){
  start = dat6[treeID, 'start']
  end = dat6[treeID, 'end']
  snp = start+1
  
  cat(treeID,'\n')
  
  filename_plot = paste('./treeViews/chr6_',geneName,'/chr6_trID',treeID,'_start',start,'_end',end,'_',geneName, sep='')
  plt = drawTree(PATH_TO_RELATE, filename_haps, filename_sample, filename_anc, filename_mut, 
                 filename_poplabels, years_per_gen, snp, filename_plot, geneName,
                 makeFiles = FALSE)
  
  rosList[[which(treeIDs == treeID)]] = plt
}

par(mfrow=c(1,1))
rosTr2 = read.tree(text = as.character(dat6[which(dat6$treeID == 513856), 'newickTree']))
# plot.phylo(rosTr2, direction='right', cex=0.5)
plot.phylo(rosTr2, tip.color = tipInfo[match(as.integer(rosTr2$tip.label), tipInfo$IDs), "tipColors"], cex=0.5, direction='down')

rosTr3 = read.tree(text = as.character(dat6[which(dat6$treeID == 513857), 'newickTree']))
# plot.phylo(rosTr3, direction='right', cex=0.5)
plot.phylo(rosTr3, tip.color = tipInfo[match(as.integer(rosTr3$tip.label), tipInfo$IDs), "tipColors"], cex=0.5, direction='down')





## Chr2: Flavia ----
geneName = "flavia"

PATH_TO_RELATE = '~/_softwares/relate_v1.2.2'
filename_haps = './Chr2/Chr2.haps.gz'
filename_sample = './Chr2/Chr2.sample'
filename_poplabels = './AvePla.MY.n74.poplabels'
filename_anc = './Chr2/rel_chr2_joint.anc.gz'
filename_mut = './Chr2/rel_chr2_joint.mut.gz'
years_per_gen = 3
snpPos = c(53652042, 53652180, 53710256, 53710327, 53712795, 53712984, 53714175)

flaviaList = list()
for (snp in snpPos){
    cat(snp, '\n')
    filename_plot = paste('./treeViews/chr2_',snp,'_',geneName, sep='')
    plt = drawTree(PATH_TO_RELATE, filename_haps, filename_sample, filename_anc, filename_mut,
                   filename_poplabels, years_per_gen, snp, filename_plot, geneName,
                   makeFiles = TRUE)

    flaviaList[[which(snpPos == snp)]] = plt
}


## Chr5: Rubia ----
geneName = "rubia"

PATH_TO_RELATE = '~/_softwares/relate_v1.2.2'
filename_haps = './Chr5/Chr5.haps.gz'
filename_sample = './Chr5/Chr5.sample'
filename_poplabels = './AvePla.MY.n74.poplabels'
filename_anc = './Chr5/rel_chr5_joint.anc.gz'
filename_mut = './Chr5/rel_chr5_joint.mut.gz'
years_per_gen = 3
snpPos = c(53652042, 53652180, 53710256, 53710327, 53712795, 53712984, 53714175)

flaviaList = list()
for (snp in snpPos){
  cat(snp, '\n')
  filename_plot = paste('./treeViews/chr5_',snp,'_',geneName, sep='')
  plt = drawTree(PATH_TO_RELATE, filename_haps, filename_sample, filename_anc, filename_mut,
                 filename_poplabels, years_per_gen, snp, filename_plot, geneName,
                 makeFiles = TRUE)
  
  flaviaList[[which(snpPos == snp)]] = plt
}