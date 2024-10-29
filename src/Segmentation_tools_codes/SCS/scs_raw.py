import os
import pandas as pd
import PIL
from PIL import Image 
from src import scs
import sys

# Change max dimension for pixels
PIL.Image.MAX_IMAGE_PIXELS = 999999999

# Add here your path for SCS-main GitHub repository folder
module_path = os.path.abspath('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/scs/SCS-main')
if module_path not in sys.path:
    sys.path.insert(0, module_path)

##########
# Input data
bin_file = output_path + "/" + "transcripts_full.csv"
image_file = output_path + "/" + "cell_image_full.tiff"

##########
# Model run
# Change patch_size by desired patches x and y shapes
scs.segment_cells(bin_file, image_file, align='rigid', patch_size=1200)
