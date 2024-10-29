######
#Prepare input data to the format baysor wants
######
import os
from os import listdir
from os.path import isfile, join
import pandas as pd # type: ignore
import multiprocessing as mp

# #Running Julia executable and storing segmentation outcome LINEARLY
# #Get all csv input files required for each patch

###
# Check in output folder those already done and remove from list to continue (in case run stopped unexpectedly)
# output_breast_dataset is the output folder for baysor run
# For breast dataset
mypath = "/ibex/user/iribarxm/spatial_transcriptomics/baysor/output_breast_dataset"

#Select "baysor" if all transcripts will be used else use "baysor_hvg" if you generated an alternative version of the transcripts files with tag _hvg_
tool = "baysor"

if tool == 'baysor_hvg':
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files = [mypath + "/" + s  for s in files]
    files_check = list(filter(lambda k: '_hvg_' in k, files))
    #
    files_names = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files_names_check = list(filter(lambda k: '_hvg_'  in k, files_names))
else:
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files = [mypath + "/" + s  for s in files]
    files_check = list(filter(lambda k: '_bidcell_' not in k, files))
    #
    files_names = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files_names_check = list(filter(lambda k: '_bidcell_' not in k, files_names))

# Set folder where all input patched transcripts files are present
# For breast dataset
mypath = "/ibex/user/iribarxm/spatial_transcriptomics/datasets/breast_data/crop_transcripts"

if tool == 'baysor_hvg':
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files = [mypath + "/" + s  for s in files]
    files = list(filter(lambda k: '_hvg_'  in k, files))
    #
    files_names = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files_names = list(filter(lambda k: '_hvg_'  in k, files_names))
else:
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files = [mypath + "/" + s  for s in files]
    files = list(filter(lambda k: '_bidcell_' not in k, files))
    #
    files_names = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files_names = list(filter(lambda k: '_bidcell_' not in k, files_names))


#Remove those already analyzed (in case run stopped unexpectedly)
a = [s.split("cell_transcripts_")[1] if "cell_transcripts_" in s else s for s in files]
a = list(map(lambda i: i[:-4], a))
c = [s.split("output/")[1] if "output/" in s else s for s in files_check]
c = list(map(lambda i: i[:-4], c))
a_remaining = [item for item in a if item not in c]

def get_matching_paths(full_paths, partial_names):
    # Dictionary to map partial names to full paths
    path_dict = {}
    #
    # Populate the dictionary
    for path in full_paths:
        # Extract the partial name from the full path
        partial_name = path.split('/')[-1].replace('cell_transcripts_', '').replace('.csv', '')
        # Check for duplicates
        if partial_name in path_dict:
            print(f"Duplicate found for partial name: {partial_name}")
        else:
            path_dict[partial_name] = path
    #
    # List to store matching paths
    matching_paths = []
    #
    # Get paths that match the partial names
    for name in partial_names:
        if name in path_dict:
            matching_paths.append(path_dict[name])
        else:
            print(f"Partial name not found in full paths: {name}")
    #
    return matching_paths

filtered_list = get_matching_paths(files, a_remaining)

files = filtered_list

if tool == 'baysor_hvg':
    files_names = [s.split("cell_transcripts_hvg_")[1] if "cell_transcripts_hvg_" in s else s for s in files]
else:
    files_names = [s.split("cell_transcripts_")[1] if "cell_transcripts_" in s else s for s in files]   

# Run Baysor for each of the input patches (in parallel, n_jobs > 1 or linearly n_jobs = 1)
# For breast dataset
def multy_run(files,files_names):
    print("RUNNING: ", files_names[:-4]) #OLD [17:-4]
    os.system(f'mkdir /ibex/user/iribarxm/spatial_transcriptomics/baysor/output_temp_breast_dataset/{files_names[:-4]}') #OLD [17:-4]
    os.system(f'~/.julia/bin/baysor run -p -c /ibex/user/iribarxm/spatial_transcriptomics/baysor/data/starmap.toml -o /ibex/user/iribarxm/spatial_transcriptomics/baysor/output_temp_breast_dataset/{files_names[:-4]} -p {files}') #OLD [17:-4] -n 20
    os.system(f'cp /ibex/user/iribarxm/spatial_transcriptomics/baysor/output_temp_breast_dataset/{files_names[:-4]}/segmentation.csv /ibex/user/iribarxm/spatial_transcriptomics/baysor/output_breast_dataset/{files_names[:-4]}".csv"') #OLD [17:-4] [17:-4]

list1 = files
list2 = files_names

input = []
for n in range(len(list1)):
    inner_array = []
    inner_array.append([list1[n], list2[n]])
    input.extend(inner_array)

n_jobs=1
pool = mp.Pool(processes=n_jobs)
pool.starmap(multy_run, input) 
pool.close()
pool.terminate()
del(pool) 
