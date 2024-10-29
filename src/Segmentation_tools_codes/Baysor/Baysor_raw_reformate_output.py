import os
from os import listdir
from os.path import isfile, join
import numpy as np 
import scanpy as sc 
import anndata as ad 
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy.stats import pearsonr 
from scipy.stats import kruskal
from scipy.sparse import lil_matrix, csr_matrix, vstack 
from sys import argv
import re
import pickle

# Set folder containing segmentation outcome files for the partches analyzed
mypath_2 = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/baysor/output_breast_dataset"

files = [f for f in listdir(mypath_2) if isfile(join(mypath_2, f))]
files = [mypath_2 + "/" + s  for s in files]
files = list(filter(lambda x : re.search(r".csv",x) ,files))

# This bucle will give outcome files the format needed for unified analysis
cell_seg2_all = files
for i in range(len(cell_seg2_all)):
    sub1 = "output_breast_dataset/"
    sub2 = ":"
    a = ''.join(cell_seg2_all[i].split(sub1,2)[1].split(sub2)[0])
    sub1 = ":"
    sub2 = "_"
    b = ''.join(cell_seg2_all[i].split(sub1,1)[1].split(sub2)[0])
    c = ''.join(cell_seg2_all[i].split(sub1,2)[1].split(sub2)[1])
    #
    cell_seg2 = pd.read_csv(cell_seg2_all[i])
    cell_seg2 = cell_seg2.dropna()
    #
    cell_seg2
    cell_seg2['x_y'] = cell_seg2["x"].astype(float).astype(int).astype(str) + ':' + cell_seg2["y"].astype(float).astype(int).astype(str)
    cell_seg2["cell"] = cell_seg2["cell"].str.split('-', n=1).str[1]
    cell_seg2['Weight'] = cell_seg2["cell"].astype(int)
    cell_seg2 = cell_seg2[["x_y","Weight"]]
    #
    # Set a folder to store reformated segmentation files
    cell_seg2.to_csv(os.path.join('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/baysor/output_breast_dataset/reformated_output_breast/'+ str(a) + ':' + str(b) + '_' + str(c) + ':' + str(c) + '.txt'), index=False, sep='\t', header=None)


