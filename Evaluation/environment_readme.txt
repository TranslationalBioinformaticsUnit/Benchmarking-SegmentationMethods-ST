#create the environment
conda env create -f /ibex/user/iribarxm/general_kernel_versions/general_kernel.yml

#activate this environment before running cellpose
conda activate kernel

#If needed to remove
conda env remove -n kernel

#Later install bidcell
#python -m pip install kernel