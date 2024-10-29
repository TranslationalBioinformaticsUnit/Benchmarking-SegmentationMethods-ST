import numpy as np 
import matplotlib.pyplot as plt 
from skimage import io, color, filters 
import spateo as st 
import scanpy as sc 
import os
from os import listdir
from os.path import isfile, join
import re

# Run segmentation for a image
def process_image(image_path):
    # Load the image
    image = io.imread(image_path)
    #
    # If the image is in color, convert it to grayscale
    if len(image.shape) == 3:
        image_gray = color.rgb2gray(image)
    else:
        image_gray = image
    #
    # Convert the image to an AnnData object
    adata = sc.AnnData(X=image_gray)
    #
    # Set the required __type key in the .uns attribute
    adata.uns['__type'] = 'AGG'
    #
    # Add the image data to a layer named 'stain'
    adata.layers['stain'] = image_gray
    #
    # Nucleus segmentation from staining image
    st.cs.mask_nuclei_from_stain(adata, otsu_classes=4, otsu_index=1)
    st.cs.find_peaks_from_mask(adata, 'stain', 7)
    #
    # Run the watershed algorithm
    st.cs.watershed(adata)
    #
    return adata

#Pipeline for watershed segmentation outcome generation
def postprocess(startx, starty, patchsize, image_path):
    #Generate segmentation results
    adata = process_image(image_path)
    #
    startx = str(startx)
    starty = str(starty)
    patchsize = str(patchsize)
    #
    #Create empty file first
    file_path = f'{myoutputpath}/watershed_result/spot2nucl_{startx}:{starty}:{patchsize}:{patchsize}.txt'
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, 'w') as fw:
    # The file is created and empty
        pass
    #
    #Export nucleus watershed segmentation
    nucl_labels = adata.layers['stain_labels']
    fw = open(myoutputpath + '/watershed_result/spot2nucl_' + startx + ':' + starty + ':' + patchsize + ':' + patchsize + '.txt', 'w')
    for i in range(nucl_labels.shape[0]):
        for j in range(nucl_labels.shape[1]):
            if nucl_labels[i, j] > 0:
                fw.write(str(i) + ':' + str(j) + '\t' + str(nucl_labels[i, j]) + '\n')
    fw.close()

# Generate outputs for all patches
# Load the image

#SET PATH FOR IMAGES TO BE SEGMENTED USING WATERSHED
mypath = "/ibex/user/iribarxm/spatial_transcriptomics/datasets/breast_data/crop_images"
myoutputpath = "/ibex/user/iribarxm/spatial_transcriptomics/watershed"

files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
files = [mypath + "/" + s  for s in files]
files = list(filter(lambda k: '.png' in k, files))

for file in files:
    #Getting values from file name
    sample_n = re.sub(".*/", "", file)
    sample_n = re.sub(".*image_(.*?)\\..*", "\\1", sample_n)
    x_ind = re.sub(r":.*", "", sample_n)
    y_ind = re.sub(".*:(.*?)_.*", "\\1", sample_n)
    patch_size_ind = re.sub(r".*:(.*?)(?:_|$).*", r"\1", sample_n)
    print("Running Watershed in sample: " +  sample_n)
    #
    image_path = file
    postprocess(x_ind,y_ind,patch_size_ind, image_path)