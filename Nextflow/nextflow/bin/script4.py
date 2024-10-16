import os
import pandas as pd
import PIL
from PIL import Image
import sys

# Adjust the module path to include the SCS package directory
module_path = os.path.abspath('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/scs/SCS-main')
if module_path not in sys.path:
    sys.path.insert(0, module_path)

from src import scs  # Import the scs module for segmentation

# This is needed for very large image files
PIL.Image.MAX_IMAGE_PIXELS = 999999999

# Set up the inputs from the Nextflow pipeline
bin_file = sys.argv[1]  # This will be the 'transcripts_formated.csv' file
image_file = sys.argv[2]  # This will be the 'Image_formated.tiff' file
output_path = sys.argv[3]  # Output path to store results

# Ensure the output directory exists
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Log the input parameters
print(f"Running segmentation with the following inputs:")
print(f"Transcript file: {bin_file}")
print(f"Image file: {image_file}")
print(f"Output path: {output_path}")

# Run the SCS cell segmentation model
scs.segment_cells(
    bin_file, 
    image_file, 
    align='rigid',  # Specify alignment type ('rigid')
    patch_size=1200  # Set patch size (this can be adjusted)
)
