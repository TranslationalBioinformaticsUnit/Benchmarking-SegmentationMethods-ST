import os
import pandas as pd  # type: ignore
import PIL  # type: ignore
from PIL import Image  # type: ignore
import sys

# Define the path for the SCS module
module_path = os.path.abspath('/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/scs/SCS-main')
if module_path not in sys.path:
    sys.path.insert(0, module_path)

from src import scs  # Ensure your SCS module is properly imported

PIL.Image.MAX_IMAGE_PIXELS = 999999999

# Read input parameters from command line
transcripts_input = sys.argv[1]
image_input = sys.argv[2]
output_path = sys.argv[3]
color_code = sys.argv[4]

# Load the transcripts matrix from CSV
transcripts_matrix = pd.read_csv(os.path.join(transcripts_input), sep=',')
transcripts_matrix = transcripts_matrix.rename(columns={'target': 'geneID', 'z_raw': 'MIDCounts', 'cell_ID': 'cell_id'})
transcripts_matrix['MIDCounts'] = 1

# Scale x,y transcripts axis to image pixel axis
needed = "true"  # Adjust as necessary

if needed == "true":
    # Subset columns
    cols = ['geneID', 'x', 'y', 'MIDCounts', 'cell_id']
    transcripts_matrix = transcripts_matrix[cols]
    #
    # Load the dataset corresponding image
    im = Image.open(image_input)
    # Trasform image to gray scale
    if color_code == 'rgb':
        im = im.convert("L")
    imgwidth, imgheight = im.size
    #
    def rescale_vector(vec, new_min, new_max):
        old_min = min(vec)
        old_max = max(vec)
        #
        # Rescale the vector
        new_vec = [(x - old_min) / (old_max - old_min) * (new_max - new_min) + new_min for x in vec]
        return new_vec
    #
    # Replace them
    rescaled_vec_x = rescale_vector(transcripts_matrix["x"], 1, imgwidth)
    rescaled_vec_y = rescale_vector(transcripts_matrix["y"], 1, imgheight)
    #
    transcripts_matrix["x"] = rescaled_vec_x
    transcripts_matrix["y"] = rescaled_vec_y
    #
    # Export the file
    transcripts_matrix["x"] = transcripts_matrix["x"].astype(int)
    transcripts_matrix["y"] = transcripts_matrix["y"].astype(int)
    transcripts_matrix["MIDCounts"] = transcripts_matrix["MIDCounts"].astype(int)
    #
    os.makedirs(output_path, exist_ok=True)
    #
    # Export full with cellID
    transcripts_matrix.to_csv(os.path.join(output_path, "transcripts_formated_full.csv"), index=False, sep=',', header=True)
    #
    # Export partial with no cellID
    cols = ['geneID', 'x', 'y', 'MIDCounts']
    transcripts_matrix = transcripts_matrix[cols]
    transcripts_matrix.to_csv(os.path.join(output_path, "transcripts_formated.csv"), index=False, sep='\t', header=True)
    #
    # Export image in TIFF format
    im.save(os.path.join(output_path, 'Image_formated.tiff'), format='TIFF')