##### Python script for tsinfer & tsdate
##### author: Arka Pal
##### written: 22.07.2024
##### update: 22.07.2024

##### USAGE: python run_ts.py <options>

##### $1 - VCF file
##### $2 - ancestral file
##### $3 - sample_data_file
##### $4 - tsinfer trees
##### $5 - tsdate trees
##### $6 - tsdate trees


## Import modules
## -----
print("Importing modules")
import tsinfer
import tsdate
import tskit
import cyvcf2
import numpy as np
import os
import pandas as pd
import msprime
from IPython.display import SVG, display
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys
from datetime import datetime
print("Done----\n")


## Read input
## -----
print("Reading input")
vcf_path = sys.argv[1]
ancestral = sys.argv[2]
sample_data_file = sys.argv[3]
tsinfer_tree_file = sys.argv[4]
tsdate_tree_file = sys.argv[5]
treeList_file = sys.argv[6]
print("Done----\n")


## Set directory
## -----
os.chdir('/nfs/scistore18/bartogrp/apal/snap_hap_repHZ/tsinfer')


## Function: Define function to add sites into samples for tsinfer
## -----
print("Define functions")
def add_diploid_sites(vcf, ancestral, samples, L):

    """
    Read ancestral alleles file and extract anecstral alleles
    """
    allelePolarised = pd.read_csv(ancestral)
    allelePolarised = allelePolarised.loc[:,['pos', 'chrom', 'ref', 'alt', 'anc_allele', 'der_allele']]
    allelePolarised['anc'] = (allelePolarised.loc[:,'ref'] != allelePolarised.loc[:,'anc_allele']).astype(int)
    
    """
    Read the sites in the vcf and add them to the samples object.
    """
    allele_chars = set("ATGCatgc*")
    pos = 0
    siteID = 0

    # bar = progressbar.ProgressBar().start()
    progressBar = tqdm(total = L, desc="Reading sites", unit="iter")
    for variant in vcf:  # Loop over variants, each assumed at a unique site    
        # bar.update(siteID)
        progressBar.update(1)
        allele_chars = set("ATGCatgc*")
        
        pos = variant.POS
        
        alleles = [variant.REF.upper()] + [v.upper() for v in variant.ALT]
        allele_ancestral = allelePolarised.loc[allelePolarised['pos'] == pos, 'anc'].values[0] if pos in allelePolarised['pos'].values else 0
        
        # Check we have ATCG alleles
        for a in alleles:
            if len(set(a) - allele_chars) > 0:
                print(f"Ignoring site at pos {pos}: allele {a} not in {allele_chars}")
                continue
        
        # Map original allele indexes to their indexes in the new alleles list.
        genotypes = [g for row in variant.genotypes for g in row[0:2]]

        samples.add_site(pos, genotypes, alleles, ancestral_allele=allele_ancestral)
        siteID += 1
    # bar.finish()
    progressBar.close()
print("Done----\n")


## Variables
## -----
demes=['AveFR','AveY', 'PlaFR', 'PlaY']
demeSize = [19, 19, 18, 18]
popList = np.repeat(demes, demeSize)
vcf = cyvcf2.VCF(vcf_path)
L = sum([1 for _ in vcf])
print('Total no. of sites: {}\n'.format(L))


## Read vcf and load into tsinfer
## -----
print('Read and load vcf')
print('vcf:{}'.format(vcf_path))
print('ancestral allele data:{}'.format(ancestral))
start_time = datetime.now()
vcf = cyvcf2.VCF(vcf_path)
with tsinfer.SampleData(path = sample_data_file) as sample_data:
    ## Define populations
    sample_data.add_population(metadata={"name": "AveFR"})
    sample_data.add_population(metadata={"name": "PlaY"})
    sample_data.add_population(metadata={"name": "PlaFR"})
    sample_data.add_population(metadata={"name": "PlaY"})

    ## Define inidividuals
    for sampleName, pop in zip(vcf.samples, popList):
        popIndex = demes.index(pop)
        sample_data.add_individual(ploidy=2, population=popIndex, metadata={"names":sampleName})

    ## Add sites and genotypes
    add_diploid_sites(vcf, ancestral, sample_data, L)
print(
    "Sample file created for {} samples ".format(sample_data.num_samples)
    + "({} individuals) ".format(sample_data.num_individuals)
    + "with {} variable sites.".format(sample_data.num_sites),
    flush=True,
)
print(f"Done in {datetime.now()-start_time} seconds ----\n")


## Infer trees
## -----
print("Infer trees")
start_time = datetime.now()
# sample_data = tsinfer.load(sample_data_file)

ts_AvePlaFrYe = tsinfer.infer(sample_data, 
                              num_threads=8, 
                              path_compression=True, 
                              progress_monitor=True)
ts_AvePlaFrYe.dump(tsinfer_tree_file)
# ts_AvePlaFrYe = tskit.load(tsinfer_tree_file)
print(f"tsinfer object saved at {tsinfer_tree_file}")
print(f"Done in {datetime.now()-start_time} seconds ----\n")


## Simplify trees
print("Remove unary nodes")
start_time = datetime.now()
tsSimp_AvePlaFrYe = ts_AvePlaFrYe.simplify(keep_unary=False)
print(f"Done in {datetime.now()-start_time} seconds ----\n")


## Date trees
print("tsdate")
start_time = datetime.now()
mu = 5.7e-9
Ne = tsSimp_AvePlaFrYe.diversity()/(4*mu)
tsD_AvePlaFrYe = tsdate.date(tsSimp_AvePlaFrYe, mutation_rate=mu,
                             method='variational_gamma')
tsD_AvePlaFrYe.dump(tsdate_tree_file)
print(f"Done in {datetime.now()-start_time} seconds ----\n")


## Make newick tree list
print("Make newick tree list")
start_time = datetime.now()
treeStart = []
treeEnd = []
treeSpan = []
treeList = []
n = ts_AvePlaFrYe.num_trees
print(f"No. of trees: {n}")
progressBar = tqdm(total = n, desc="Reading trees", unit="trees")
for tree in ts_AvePlaFrYe.trees():
    progressBar.update(1)
    if (tree.interval.left > 0):
        treeStart.append(tree.interval.left)
        treeEnd.append(tree.interval.right)
        treeSpan.append(tree.span)
        # treeList.append(tree.as_newick(root=tree.roots))
        treeList.append(tree.as_newick(root=tree.root))
progressBar.close()

tsDF = pd.DataFrame({'treeStart': treeStart,
                      'treeEnd': treeEnd,
                     'treeSpan': treeSpan,
                     'tree': treeList})
tsDF.to_csv(treeList_file, sep="\t", quoting = None, index = False, encoding = 'utf-8')
print(f"Done in {datetime.now()-start_time} seconds ----\n")