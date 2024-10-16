import os
import sys
import PIL
from PIL import Image
import math
import numpy as np
import tifffile as tiff
import csv  # Import csv module if needed

import logging
logging.basicConfig(level=logging.INFO)

PIL.Image.MAX_IMAGE_PIXELS = 999999999

# Read input parameters from the command line
bin_file = sys.argv[1]
input_image = sys.argv[2]
patch_size = int(sys.argv[3])
output_path = sys.argv[4]
subset = sys.argv[5]
color_code = sys.argv[6]

#Lets get x min/max and y min/max pixels to reshape the image
def detect_min_max_x_y(bin_file,patch_size):
    r_all = [] 
    c_all = []
    with open(bin_file) as fr:
        header = fr.readline()
        for line in fr:
            _, r, c, _ = line.split()
            r_all.append(int(r))
            c_all.append(int(c))
    rmax = np.max(r_all) - np.min(r_all)
    cmax = np.max(c_all) - np.min(c_all)
    n_patches = math.ceil(rmax / patch_size) * math.ceil(cmax / patch_size)
    print(str(n_patches) + ' patches will be processed.')
    output= [rmax, np.min(r_all), cmax, np.min(c_all)]
    return output

#Lets crop into n patches os size 1200x1200 base on previusly x y detected
def crop(path, input, height, width, min_max_x_y, subset, color_code):
    # Create the output directory if it doesn't exist
    os.makedirs(path, exist_ok=True)
        
    im = Image.open(input)
    #Trasform image to gray scale
    if color_code == 'rgb':
        im = im.convert("L")
    #
    if subset == "True":
        #Subset only the interest part of image
        box = (min_max_x_y[3], min_max_x_y[1], min_max_x_y[2]+min_max_x_y[3], min_max_x_y[0]+min_max_x_y[1])
        im = im.crop(box)
    # Create supplementary directory for full images
    full_image_path = os.path.join(path, 'full_image')
    os.makedirs(full_image_path, exist_ok=True)
    im.save(os.path.join(full_image_path,"cell_image_full.png"))
    im.save(os.path.join(full_image_path,"cell_image_full.tif"))
    #    
    imgwidth, imgheight = im.size
    for i in range(0,imgheight,height):
        for j in range(0,imgwidth,width):
            print("Patch_section"+" "+str(i)+" "+str(j))
            box = (j, i, j+width, i+height)
            a = im.crop(box)
            #Export
            a.save(os.path.join(path,"cell_image_"+str(i)+":"+str(j)+"_"+str(width)+":"+str(height)+".png"))
            a.save(os.path.join(path,"cell_image_bidcell_"+str(i)+":"+str(j)+"_"+str(width)+":"+str(height)+".tif"))

# Run the detection of min/max x/y
min_max_x_y = detect_min_max_x_y(bin_file, patch_size)

# Crop the image based on the parameters provided
crop(output_path, input_image, patch_size, patch_size, min_max_x_y, subset, color_code)
