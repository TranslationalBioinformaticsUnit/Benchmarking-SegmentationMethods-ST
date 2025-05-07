#create the environment
conda env create -f /ibex/user/iribarxm/spatial_transcriptomics/scs/environment_scs.yml

#activate this environment before running scs
conda activate scs

#If needed to remove
conda env remove -n scs