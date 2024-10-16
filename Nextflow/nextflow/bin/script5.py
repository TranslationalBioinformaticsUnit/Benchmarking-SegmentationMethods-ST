import numpy as np
import matplotlib.pyplot as plt
from cellpose import denoise, models, io
from cellpose import plot
from cellpose.io import imread
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import sys

# Logging setup for cellpose
io.logger_setup()

# Command-line arguments from Nextflow
patches_dir = sys.argv[1]  # Directory containing image patches
cellpose_model = sys.argv[2]
cellpose_restore_type = sys.argv[3]
output_path = sys.argv[4]

# Ensure the output directory exists
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Select the pre-trained model to use
model = denoise.CellposeDenoiseModel(gpu=True, model_type=cellpose_model, restore_type=cellpose_restore_type)

# List of image patch files in the specified directory
files = [f for f in listdir(join(patches_dir, output_path)) 
         if isfile(join(patches_dir, output_path, f)) and '_bidcell_' not in f]
files = [join(patches_dir,output_path, f) for f in files]

# Load the images from the file paths
imgs = [imread(f) for f in files]
nimg = len(imgs)

# Define CHANNELS for segmentation
channels = [[0, 0]]  # Assuming grayscale images

# Perform segmentation on the image patches
masks, flows, styles, diams = model.eval(imgs, diameter=None, channels=channels)

# Process and save the results
for idx in range(nimg):
    maski = masks[idx]
    flowi = flows[idx][0]
    
    # Plot and save the segmentation results as a PDF
    fig = plt.figure(figsize=(12, 5))
    plot.show_segmentation(fig, imgs[idx], maski, flowi, channels=channels[0])
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "cellpose_spot2cell_" + files[idx].split('/')[-1][10:-4] + ".pdf"), format="pdf", bbox_inches="tight")
    
    # Generate a spot2cell matrix from the segmentation mask
    roi_table = pd.DataFrame(masks[idx]).stack().reset_index().rename(columns={'level_0': 'Source', 'level_1': 'Target', 0: 'Weight'})
    roi_table['x_y'] = roi_table['Source'].astype(str) + ':' + roi_table['Target'].astype(str)
    roi_table['Weight'] = roi_table['Weight'].astype(float)
    
    # Subset and reorder columns, and remove pixels not assigned to any cell
    roi_table = roi_table[['x_y', 'Weight']]
    roi_table = roi_table[roi_table['Weight'] != 0.0]
    
    # Save the spot2cell matrix as a TXT file
    roi_table.to_csv(os.path.join(output_path, "cellpose_spot2cell_" + files[idx].split('/')[-1][10:-4] + ".txt"), index=False, sep='\t', header=None)
