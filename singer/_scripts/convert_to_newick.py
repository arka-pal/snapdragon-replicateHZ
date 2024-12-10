#!/usr/bin/env python3

##### USAGE: module load python; python convert_to_newick.py <args>
##### $1: treeSeqFile
##### $2: outNewickFile


## import libraries
import tskit
import pandas as pd
import sys
from tqdm import tqdm

## Read arguments
# # Check if the correct number of arguments is provided
# if len(sys.argv) != 3:
#     print("Usage: ./script.py <var1> <var2>")
#     sys.exit(1)
treeSeqFile = sys.argv[1]
outNewickFile = sys.argv[2]

## Print argurments
print('Tree Sequence File:', treeSeqFile)
print('output File:', outNewickFile)

## Read tree sequence
tr = tskit.load(treeSeqFile)
# tr

## Convert tree sequnece to newick format
n = tr.num_trees
progressBar = tqdm(total = n, desc="Reading trees", unit="trees")
treeStart = []
treeEnd = []
treeSpan = []
treeList = []
for tree in tr.trees():
    progressBar.update(1)
    if (tree.interval.left > 0):
        treeStart.append(tree.interval.left)
        treeEnd.append(tree.interval.right)
        treeSpan.append(tree.span)
        # treeList.append(tree.as_newick(root=tree.roots))
        treeList.append(tree.as_newick(root=tree.root))
progressBar.close()

## Make dataframe with newick format trees
trDF = pd.DataFrame({'treeStart': treeStart,
                      'treeEnd': treeEnd,
                     'treeSpan': treeSpan,
                     'tree': treeList})

## Write File
trDF.to_csv(outNewickFile, sep="\t", quoting = None, index = False, encoding = 'utf-8')