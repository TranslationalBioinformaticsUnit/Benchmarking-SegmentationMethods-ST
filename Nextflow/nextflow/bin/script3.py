import os
import sys
import pandas as pd

# Read input parameters from the command line
bin_file = sys.argv[1]
output_path = sys.argv[2]
patch_size = int(sys.argv[3])

# Create the output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)

def main(bin_file, output_path):
    # Load RNA matrix
    rna_matrix = pd.read_csv(bin_file, sep=',')
    rna_matrix["MIDCounts"] = 1
    rna_matrix['z'] = 0
    #
    # List of tools to process
    tools = ['baysor', 'BIDCell']
    #
    for tool in tools:
        print(f"Processing tool: {tool}")
        #
        if tool == 'baysor' or tool == 'baysor_hvg':
            # For baysor
            rna_matrix_tool = rna_matrix.iloc[:, [1, 2, 0, 5]]
            rna_matrix_tool.columns = ['x', 'y', 'gene', 'z']
            #
        if tool == 'BIDCell':
            # For BIDCell
            rna_matrix_tool = rna_matrix.iloc[:, [1, 2, 0, 3, 5]]
            rna_matrix_tool.columns = ['x', 'y', 'geneID', 'MIDCounts', 'z']
        #
        # Move coordinates to (0,0)
        rna_matrix_tool['x'] -= rna_matrix_tool['x'].min()
        rna_matrix_tool['y'] -= rna_matrix_tool['y'].min()
        #
        # Subset into patches
        imgheight = rna_matrix_tool['y'].max()
        imgwidth = rna_matrix_tool['x'].max()
        height = patch_size
        width = patch_size
        count = 0
        #
        for y in range(0, int(imgheight), height):
            for x in range(0, int(imgwidth), width):
                count += 1
                print(f"Processing patch {count} for {tool}")
                #
                # Subset rows within the current patch
                temp_y = rna_matrix_tool[rna_matrix_tool['y'].between(y, y + height - 1)]
                #
                if temp_y['x'].min() < (x + width):
                    temp_y = temp_y[temp_y['x'].between(x, x + width - 1)]
                    #
                    # Move coordinates to 1200x1200 range
                    temp_y['y'] -= y
                    temp_y['x'] -= x
                    #
                    if not temp_y.empty: 
                        # Export the current patch based on the tool
                        if tool == 'BIDCell':
                            temp_y.to_csv(os.path.join(output_path, f'cell_transcripts_bidcell_{y}:{x}_{height}:{width}.csv'), index=False, sep=',', header=True)
                        elif tool == 'baysor':
                            temp_y.to_csv(os.path.join(output_path, f'cell_transcripts_{y}:{x}_{height}:{width}.csv'), index=False, sep=',', header=True)
                        elif tool == 'baysor_hvg':
                            temp_y.to_csv(os.path.join(output_path, f'cell_transcripts_hvg_{y}:{x}_{height}:{width}.csv'), index=False, sep=',', header=True)
                    else:
                        print(f"Patch {y}:{y+height} _ {x}:{x+width} contains no cells.")

main(bin_file, output_path)