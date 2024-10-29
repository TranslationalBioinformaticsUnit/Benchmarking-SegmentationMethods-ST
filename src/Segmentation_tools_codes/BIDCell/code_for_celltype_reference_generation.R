# USE conda activate kernel (using multipropose installation of python and R)

# Required for R version 4.3.3 to install Seurat4 in order
#install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz', repos = NULL)
#install.packages('https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-60.tar.gz', repos = NULL)
#conda install conda-forge::r-igraph
#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

# Load required modules from linux (if not available in your system install it)
#module load imagemagick

library("Seurat")
library("here")
library("ggplot2")
library("dplyr")
library("reticulate")
library("anndata")
library("arrow")
use_python("/home/iribarxm/.virtualenvs/r-reticulate/bin/python", required = TRUE)
library("biomaRt")
library("Matrix.utils")
library("Matrix")
library("dplyr")

# 1st Load x,y count data from spatial transcriptomics data 
# 2nd Load it into R
# 3rd Transform it into gene x cell matrix
# 4th Annotate base on cell_type annotated data 
# 5th Cell annotation through seurat4 pipeline
# 6th Do the transformation to generate the matrix required by BIDcell and the 2 extra inputs from it

# Set working directory
setwd("/ibex/user/iribarxm/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data")

# Load data
df = read.delim("/ibex/user/iribarxm/spatial_transcriptomics/datasets/breast_data/breast_transcripts_formated_full.csv",sep=",")

# Transform it into gene x cell matrix
df$geneID = as.factor(df$geneID)
df$cell_id = as.factor(df$cell_id)
pivot_df = dMcast(df, geneID~cell_id, value.var = 'MIDCounts')

# Cell annotation through seurat4 pipeline
brain = CreateSeuratObject(pivot_df, project = "SeuratProject", assay = "RNA",
                            min.cells = 0, min.features = 0, names.field = 1,
                            names.delim = "_", meta.data = NULL)

brain = SCTransform(brain, assay = "RNA", ncells = 3000, verbose = FALSE)
brain = RunPCA(brain, assay = "SCT", verbose = FALSE)
brain = FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain = FindClusters(brain, verbose = FALSE)
brain = RunUMAP(brain, reduction = "pca", dims = 1:30)

save.image("annotation_worflow_from_ref.Rdata")

# Load cell references and cell annotation metadata
# Use this as reference data to transfer cell annotations
adata <- read_h5ad("/ibex/user/iribarxm/spatial_transcriptomics/datasets/breast_data/BC_atlas_xe.h5ad")
metadata = data.frame(adata$obs)
rownames(metadata) = gsub("_", "-", rownames(metadata))
expression = t(data.frame(adata$X))
rownames(expression) = gsub("_", "-", rownames(expression))

allen_reference_annotation = metadata
allen_reference = expression

allen_reference = as.matrix(allen_reference)
allen_reference = Matrix(allen_reference, sparse = TRUE)

allen = CreateSeuratObject(allen_reference, project = "SeuratProject", assay = "RNA",
                            min.cells = 0, min.features = 0, names.field = 1,
                            names.delim = "_", meta.data = allen_reference_annotation)

# Depends if data is already normalized "no" or not "yes"
to_normalize = "no"

if (to_normalize == "yes"){
    allen = SCTransform(allen, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
} else if (to_normalize == "no"){
    allen@assays$RNA$data = allen@assays$RNA$counts
    allen = FindVariableFeatures(allen, nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)
}

# Transfer labels
DefaultAssay(brain) = "RNA"
DefaultAssay(allen) = "RNA"

# Needed to have files in logNormalize assay
brain@assays$RNA$data = brain@assays$SCT$data
brain@assays$RNA$scale.data = brain@assays$SCT$scale.data

# Find anchors for label transfer
anchors = FindTransferAnchors(reference = allen, query = brain, normalization.method = "LogNormalize")
predictions.assay = TransferData(anchorset = anchors, refdata = allen$celltype_major, dims = 1:30)
brain = AddMetaData(brain, metadata = predictions.assay)

# Subset relevant cell_type assignation
Idents(brain) <- "predicted.id" 

# Export csv required for BIDcell and to compute the other 2 files required
DefaultAssay(object = brain) <- "RNA"
brain = NormalizeData(object = brain)

sc_brain_csv = t(AverageExpression(object = brain, group.by = c('predicted.id'))$RNA)
sc_brain_csv = as.data.frame(sc_brain_csv, optional=TRUE)
sc_brain_csv$cell_type = rownames(sc_brain_csv)
sc_brain_csv$atlas = "breast_ref"

rownames(sc_brain_csv) = NULL

## Add the cell_type index column ct_idx
cell_types = sort(unique(sc_brain_csv$cell_type))
ct_idx     = 0:(length(cell_types)-1)
names(ct_idx) = cell_types

sc_brain_csv$ct_idx = sapply(sc_brain_csv$cell_type, function(x) ct_idx[x])

write.csv(sc_brain_csv, "/ibex/user/iribarxm/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/sc_references/sc_breast.csv", row.names=TRUE) 


#Export 1000 hvg (needed for BIDCell)
brain = FindVariableFeatures(object = brain)
if (length(brain@assays$SCT@var.features) >= 1000){
  fp_selected_genes = brain@assays$SCT@var.features[1:1000]
}else{
  fp_selected_genes = brain@assays$SCT@var.features[1:length(brain@assays$SCT@var.features)]
}

# Remove unwanted gene names from list
fp_selected_genes <- rownames(brain)[!grepl("NegControlCodeword|BLANK-|NegControlProbe|antisense-", rownames(brain))]

# Export outcomes
write.table(fp_selected_genes, "/ibex/user/iribarxm/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/sc_references/fp_selected_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

save.image("annotation_worflow_from_ref.Rdata")


