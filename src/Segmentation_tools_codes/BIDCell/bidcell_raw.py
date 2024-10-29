# Runn BIDcell over all patches we generated and save the outcomes in a folder for corresponding patch in parallel
import os
from os import listdir
from os.path import isfile, isdir, join
# Change working directory to GitHub downloaded repository sourcecode from BIDCell
os.chdir('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main')
from bidcell import BIDCellModel
import pandas as pd 
import numpy as np 
import argparse
import pandas as pd 
import natsort 

# Working data path (with default yaml)
#For breast dataset
os.chdir("/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data")

# SPATIAL TRANSCRIPTS (this rules over images)
# For breast dataset
mypath = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/datasets/breast_data/crop_transcripts"

# Check already processed ones (in case run stopped unexpectedly)
current_output = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output"
folders = [f for f in listdir(current_output) if isdir(join(current_output, f))]
folders = ["cell_transcripts_bidcell_" + s + ".csv"  for s in folders]

# Check all files
files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
files_transcripts = [mypath + "/" + s  for s in files]
files_names_transcripts = [f for f in listdir(mypath) if isfile(join(mypath, f))]
files_names_transcripts = list(filter(lambda k: '_bidcell_' in k, files_names_transcripts))

# Select only remaining ones (in case run stopped unexpectedly)
files_names_transcripts = list(set(files_names_transcripts).symmetric_difference(set(folders)))

def main(config):
    #
    ref_df = pd.read_csv(config.fp_ref) #index_col=0
    n_genes = ref_df.shape[1] - 3
    print("Ref data shape", ref_df.shape)
    #
    cell_types = ref_df["cell_type"].tolist()
    cell_types = natsort.natsorted(list(set(cell_types)))
    print(cell_types)
    n_cell_types = len(cell_types)
    #
    ref_expr = ref_df.iloc[:, :n_genes].to_numpy()
    gene_names = ref_df.columns[:n_genes]
    #
    # Find genes with expressions in bottom 10% percentile for every ref cell type
    pct_10 = np.percentile(ref_expr, 10, axis=1, keepdims=True)
    pct_10 = np.tile(pct_10, (1, n_genes))
    low_expr_true = np.zeros(pct_10.shape)
    low_expr_true[ref_expr <= pct_10] = 1
    #
    # Find overlap for different ref samples of the same cell type
    ct_idx = ref_df["ct_idx"].to_numpy()
    low_expr_true_agg = np.zeros((n_cell_types, n_genes))
    for ct in range(n_cell_types):
        rows = np.where(ct_idx == ct)[0]
        low_expr_true_ct = low_expr_true[rows]
        low_expr_true_agg[ct, :] = np.prod(low_expr_true_ct, axis=0)
    #
    # Set overlaps to 0
    overlaps = np.sum(low_expr_true_agg, 0)
    too_many = np.where(overlaps > config.max_overlaps_neg)[0]
    low_expr_true_agg[:, too_many] = 0
    #
    df_neg = pd.DataFrame(low_expr_true_agg, index=cell_types, columns=gene_names)
    #
    # Find genes with expressions in top 90% percentile for every ref cell type
    pct_90 = np.percentile(ref_expr, 90, axis=1, keepdims=True)
    pct_90 = np.tile(pct_90, (1, n_genes))
    high_expr_true = np.zeros(pct_90.shape)
    high_expr_true[ref_expr >= pct_90] = 1
    #
    # Find overlap for different ref samples of the same cell type
    ct_idx = ref_df["ct_idx"].to_numpy()
    high_expr_true_agg = np.zeros((n_cell_types, n_genes))
    for ct in range(n_cell_types):
        rows = np.where(ct_idx == ct)[0]
        high_expr_true_ct = high_expr_true[rows]
        high_expr_true_agg[ct, :] = np.prod(high_expr_true_ct, axis=0)
    #
    # Set overlaps to 0
    overlaps = np.sum(high_expr_true_agg, 0)
    too_many = np.where(overlaps > config.max_overlaps_pos)[0]
    high_expr_true_agg[:, too_many] = 0
    #
    df_pos = pd.DataFrame(high_expr_true_agg, index=cell_types, columns=gene_names)
    #
    df_pos.to_csv(config.fp_pos, index=True)
    df_neg.to_csv(config.fp_neg, index=True)
    #
    print("Done")

# Change paths for your oun ones (also do this for .yaml configuration file)
def multy_run(files_names_transcripts):
    os.chdir("/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data")
    print("RUNNING: ", files_names_transcripts[25:-4])
    #Create directory for output from current patch model
    os.system(f'mkdir /ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/{files_names_transcripts[25:-4]}')
    #LINES TO MOD 7,8 AND 9
    line_number = [3,6,7,8]  # Whatever the line number you're trying to replace is
    replacement_line_3 = f'cpus: 16 '
    replacement_line_7 = f'  data_dir: /ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/{files_names_transcripts[25:-4]}'
    replacement_line_8 = f'  fp_dapi: /ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/datasets/breast_data/crop_images/cell_image_bidcell_{files_names_transcripts[25:-4]}.tif'
    replacement_line_9 = f'  fp_transcripts: /ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/datasets/breast_data/crop_transcripts/cell_transcripts_bidcell_{files_names_transcripts[25:-4]}.csv'
    #
    items = os.listdir(".")  # Gets all the files & directories in the folder containing the script
    #
    #Mod default yaml to adjust for patch specific files
    for file_name in items:  # For each of these files and directories,
        if file_name.lower().endswith(".yaml"):  # check if the file is a YAML. If it is:
            os.system(f'cp {file_name} /ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/{files_names_transcripts[25:-4]}/{files_names_transcripts[25:-4]}.yaml')
            file_name = '/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/'+files_names_transcripts[25:-4]+'/'+files_names_transcripts[25:-4]+'.yaml'
            with open(file_name, "r+", encoding='utf-8') as file:  # Safely open the file
                lines = file.read().splitlines()   # Read its contents
                lines[line_number[0]] = replacement_line_3  # Replace the line
                lines[line_number[1]] = replacement_line_7  # Replace the line
                lines[line_number[2]] = replacement_line_8  # Replace the line
                lines[line_number[3]] = replacement_line_9  # Replace the line
                data = '\n'.join(lines) + '\n'
                file.seek(0)
                file.truncate()
                file.write(data)  # And save the file

    ######
    # Alternative way for working with patches
    ######
    # Run model preprocess
    model = BIDCellModel('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/'+files_names_transcripts[25:-4]+'/'+files_names_transcripts[25:-4]+'.yaml') # type: ignore
    model.preprocess()
    #
    # Fix files and reload yaml on new files for training (base on patch genes)
    df = pd.read_csv('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/sc_references/sc_breast.csv')
    genes_to_use = pd.read_csv('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/'+files_names_transcripts[25:-4]+'/all_gene_names.txt', sep=" ", header=None)
    genes_to_use = genes_to_use[0].tolist()
    genes_to_use = genes_to_use + ['cell_type','atlas','ct_idx']
    df = df[genes_to_use] 
    df['cell_type'] = df['cell_type'].astype(str)
    df['atlas'] = df['atlas'].astype(str)
    df['ct_idx'] = df['ct_idx'].astype(str)
    df.to_csv('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/'+files_names_transcripts[25:-4]+'/sc_breast.csv', index=True)  
    #
    os.chdir('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/'+files_names_transcripts[25:-4])
    parser = argparse.ArgumentParser()
    parser.add_argument("--fp_ref",default='/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/'+files_names_transcripts[25:-4]+'/sc_breast.csv',type=str,help="ref data",)
    parser.add_argument("--max_overlaps_pos",default=4,type=int,help="no more than n cell types can have the same markers",)
    parser.add_argument("--max_overlaps_neg",default=15,type=int,help="no more than n cell types can have the same markers",)
    parser.add_argument("--fp_pos", default="markers_pos.csv", type=str, help="positive markers")
    parser.add_argument("--fp_neg", default="markers_neg.csv", type=str, help="negative markers")
    config = parser.parse_args()
    main(config)
    #
    # Determine number of steps for specific patch base on size and configuration sizes(48x48)
    affine = pd.read_csv('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/'+files_names_transcripts[25:-4]+'/affine.csv', sep='\t')  
    x = int(float(affine.iloc[3,1]))
    y = int(float(affine.iloc[4,1]))
    steps = round((round(x/48)*round(y/48)) *0.3)
    #
    # Recreate yaml for running over matching patch genes (after prep-filtering)
    # LINES TO MOD
    line_number = [9,10,11,73]  # Whatever the line number you're trying to replace is
    replacement_line_9 = f'  fp_ref: /ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/{files_names_transcripts[25:-4]}/sc_breast.csv'
    replacement_line_10 = f'  fp_pos_markers: /ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/{files_names_transcripts[25:-4]}/markers_pos.csv'
    replacement_line_11 = f'  fp_neg_markers: /ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/{files_names_transcripts[25:-4]}/markers_neg.csv'
    replacement_line_73 = f'  test_step: {steps}'
    #
    #
    items = os.listdir(".")  # Gets all the files & directories in the folder containing the script
    #
    # Mod default yaml to adjust for patch specific files (pre-filtered genes)
    for file_name in items:  # For each of these files and directories,
        if file_name.lower().endswith(".yaml"):  # check if the file is a YAML. If it is:
            file_name = '/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/'+files_names_transcripts[25:-4]+'/'+files_names_transcripts[25:-4]+'.yaml'
            with open(file_name, "r+", encoding='utf-8') as file:  # Safely open the file
                lines = file.read().splitlines()   # Read its contents
                lines[line_number[0]] = replacement_line_9  # Replace the line
                lines[line_number[1]] = replacement_line_10  # Replace the line
                lines[line_number[2]] = replacement_line_11  # Replace the line
                lines[line_number[3]] = replacement_line_73  # Replace the line
                data = '\n'.join(lines) + '\n'
                file.seek(0)
                file.truncate()
                file.write(data)  # And save the file    
    #
    # Run train
    model = BIDCellModel('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output/'+files_names_transcripts[25:-4]+'/'+files_names_transcripts[25:-4]+'.yaml') # type: ignore
    model.train()
    # Run predict
    model.predict()

# Run 
for n in files_names_transcripts:
    multy_run(n)
