#create the environment
conda env create -f /ibex/user/iribarxm/spatial_transcriptomics/cellpose/environment_cellpose.yml 

#activate this environment before running cellpose
conda activate cellpose

#If needed to remove
conda env remove -n cellpose