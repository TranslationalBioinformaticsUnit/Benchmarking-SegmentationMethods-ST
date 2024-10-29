import numpy as np 
import matplotlib.pyplot as plt 
from cellpose import denoise, models, io
from cellpose.io import imread 
import os
from os import listdir
from os.path import isfile, join
from cellpose import plot
import pandas as pd 

io.logger_setup()

# Set working directory
os.chdir('/ibex/user/iribarxm/spatial_transcriptomics/cellpose')

# Models avaialble: model_type='cyto' or 'nuclei' or 'cyto2' or 'cyto3'
model = denoise.CellposeDenoiseModel(gpu=True, model_type="cyto3",restore_type="denoise_cyto3")

# List of files: PUT PATH TO YOUR FILES HERE!
#For pancreas_data dataset
mypath = "/ibex/user/iribarxm/spatial_transcriptomics/datasets/pancreas_data/crop_images"

files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
files = [mypath + "/" + s  for s in files]

imgs = [imread(f) for f in files]
nimg = len(imgs)

# Define CHANNELS to run segementation on
# Grayscale=0, R=1, G=2, B=3
# Channels = [cytoplasm, nucleus]
# If NUCLEUS channel does not exist, set the second channel to 0
channels = [[0,0]]
# IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
# channels = [0,0] # IF YOU HAVE GRAYSCALE
# channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
# channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus

# If diameter is set to None, the size of the cells is estimated on a per image basis
# You can set the average cell `diameter` in pixels yourself (recommended)
# Diameter can be a list or a single number for all images
masks, flows, styles, diams = model.eval(imgs, diameter=None, channels=channels)

#OUTCOMES
output_path = "/ibex/user/iribarxm/spatial_transcriptomics/cellpose/results"
names = [f for f in listdir(mypath) if isfile(join(mypath, f))]
nimg = len(imgs)
for idx in range(nimg):
    maski = masks[idx]
    flowi = flows[idx][0]
    #
    fig = plt.figure(figsize=(12,5))
    plot.show_segmentation(fig, imgs[idx], maski, flowi, channels=channels[0]) #channels is 0 gray /channels[idx]
    plt.tight_layout()
    plt.savefig(os.path.join(output_path,"cellpose_spot2cell"+names[idx][10:-4]+".pdf"), format="pdf", bbox_inches="tight")
    #
    #spot2cell matrix from cellpose structure
    #type: ignore
    roi_table = pd.DataFrame(masks[idx]).stack().reset_index().rename(columns={'level_0':'Source','level_1':'Target', 0:'Weight'})
    roi_table['x_y'] = roi_table['Source'].astype(str) + ':' + roi_table['Target'].astype(str)
    roi_table['Weight'] = roi_table['Weight'].astype(float)
    #Subset and reorder
    roi_table = roi_table.iloc[:,2:4]
    cols = ['x_y', 'Weight']
    roi_table = roi_table[cols]
    #Remove pixels not assigned to any cell
    roi_table = roi_table[roi_table.Weight != 0.0]
    #
    roi_table.to_csv(os.path.join(output_path,"cellpose_spot2cell"+names[idx][10:-4]+".txt"), index=False, sep='\t', header=None)

