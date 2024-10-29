#create the environment
conda env create -f /ibex/user/iribarxm/spatial_transcriptomics/bidcell/environment_bidcell.yml

#activate this environment before running bidcell
conda activate bidcell

#If needed to remove
conda env remove -n bidcell

#Later install bidcell
#python -m pip install bidcell