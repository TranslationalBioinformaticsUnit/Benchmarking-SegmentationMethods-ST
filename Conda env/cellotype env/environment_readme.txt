# install pytorch and detectron2
conda install pytorch==1.9.0 torchvision==0.10.0 cudatoolkit=11.1 -c pytorch -c nvidia
python -m pip install 'git+https://github.com/facebookresearch/detectron2.git'
# (add --user if you don't have permission)

# Compile Deformable-DETR CUDA operators
git clone https://github.com/fundamentalvision/Deformable-DETR.git
cd Deformable-DETR
cd ./models/ops
module load cuda/11.8
sh ./make.sh

# clone and install the project   
pip install cellotype

# Clone the repository:
git clone https://github.com/maxpmx/CelloType.git
cd CelloType

#Then Download the model weights:
cd models
sh download.sh
cd ..

#################################################################################################

#create the environment
conda env create -f /ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/cellotype/environment_cellotype.yml

#activate this environment before running cellotype
conda activate cellotype

#If needed to remove
conda env remove -n cellotype