#Load required libraries and modules

# CellSPA installation
#library('devtools')
#install_github("SydneyBioX/CellSPA")

# Linux modules required for loading these R packages (if not available in your system, install them)
#module load geos
#module load imagemagick

#Libraries loading
library(CellSPA)
library(ggplot2)
library(ggthemes)
library(scater)
library(SingleCellExperiment)
library(SpatialExperiment)
library(rgeos)
library(dplyr)
library(pryr)
library(tidyr)
library(parallel)
library(magick)
library(viridis)
library(cluster)
library(reshape2) 
library(fmsb)
library(Matrix.utils)
library(gridExtra)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
theme_set(theme_bw())


# Paths to input files
path_to_annotation_worflow_from_ref_Rdata = "/ibex/user/iribarxm/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/annotation_worflow_from_ref.Rdata"
path_to_formated_transcripts = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/datasets/breast_data/breast_transcripts_formated_full.csv"
path_to_watershed_segmentation_outcomes = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/watershed/results_breast_data"
path_to_bidcell_segmentation_outcomes = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/bidcell/BIDCell-main/breast_cancer_data/output"
path_to_cellpose_segmentation_outcomes = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/cellpose/results_breast_dataset"
path_to_sam_segmentation_outcomes = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/sam/reformated_sam_output_breast_dataset"
path_to_sam2_segmentation_outcomes = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/sam2/reformated_sam2_output_xenium_breast_dataset"
path_to_scs_segmentation_outcomes = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/scs/path"
path_to_scs_nucleus_segmentation_outcomes = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/scs/path"
path_to_baysor_segmentation_outcomes = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/baysor/output_breast_dataset/reformated_output_breast"

# Paths to output files
path_output_files = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/valiadation/CellSPA_R_metrics/validation_outcomes"
tissue_subpath = "/breast"

#List of segmentation tools evaluated
list_of_segmentation_tools = c("watershed","bidcell","cellpose","sam","sam2","baysor") 



# STEPS

# 1st
# Load metadata

# Breast dataset
load(path_to_annotation_worflow_from_ref_Rdata)
metadata = brain@meta.data
metadata_coord = read.delim(path_to_formated_transcripts,sep=",")

# 2nd
# Load data

#Breast dataset
df = read.delim(path_to_formated_transcripts,sep=",")
df = df[,c(2,3,1,4,5)]
colnames(df) = c("x","y","geneID","MIDCounts","cell_id")

# 3rd 
# Select tools outputs from: watershed/scs_nucleus/bidcell/cellpose/scs/sam/sam2
list_of_tools = list_of_segmentation_tools
tool_used = "ST"

# 4th
# Start processing data and metadata

# Center to 0,0 x and y axis of Global image to fit the Patches (here becouse df is also loaded in the Rdata)
df['x'] = df['x'] - min(df['x'])
df['y'] = df['y'] - min(df['y'])

# Send coord information to metadata dataframe
metadata_coord$cell_id = paste0("cell_id", metadata_coord$cell_id)

# Compute centroids
centroids = metadata_coord %>%
            group_by(cell_id) %>%
            summarize(centroid_x = mean(x),
                    centroid_y = mean(y))
rownames(centroids) = centroids$cell_id

# Set same order as present in metadata
centroids = centroids[rownames(metadata),]

metadata$coord_x = centroids$centroid_x
metadata$coord_y = centroids$centroid_y

# Transform it into gene x cell matrix
df$geneID = as.factor(df$geneID)
df$label = as.factor(df$cell_id)
sce = dMcast(df, geneID~label, value.var = 'MIDCounts') 

# Fix colnames
colnames(sce) = gsub("label", "cell_id", colnames(sce))

# Subset for overlapping cells from expression matrix and metadata
sce = sce[,rownames(metadata)]

# 5th
# Load data to generate the SpatialExpreriment object
spe = SpatialExperiment::SpatialExperiment(assay = list(counts = sce),
                                            colData =  DataFrame(metadata),
                                            #rowData = rowData(sce),
                                            spatialCoordsNames = c("coord_x","coord_y"))

.initialise_CellSPA_list = function(spe) {
    spe@metadata$CellSPA = list()
    spe@metadata$CellSPA$metrics = vector("list", 4)
    names(spe@metadata$CellSPA$metrics) = c("dataset_level",
                                             "gene_level",
                                             "cell_level",
                                             "bin_level")
    return(spe)
}
.add_dataset_metrics = function(spe) {
    spe@metadata$CellSPA$dataset_level_metrics = list(num_cells = ncol(spe),
                                                       num_genes = nrow(spe))
    spe = .add_metrics(spe, c("num_cells", "num_genes"), "dataset_level")
    return(spe)
}
.add_metrics = function(spe, metrics_name, type) {
    spe@metadata$CellSPA$metrics[[type]] = sort(unique(c(spe@metadata$CellSPA$metrics[[type]],
                                                          metrics_name)))
    return(spe)
}
.cal_baseline = function(spe, metrics = c("total_transcripts", "total_genes",
                                           "total_cells", "meanExprsPct_cells")) {
    if (!"CellSPA" %in% names(spe@metadata)) {
        spe@metadata$CellSPA = list()
    }
    if ("total_transcripts" %in% metrics) {

        total_transcripts = Matrix::colSums(counts(spe))
        colData(spe)$total_transcripts = total_transcripts
        spe@metadata$CellSPA$metrics = c(spe@metadata$CellSPA$metrics,
                                          "total_transcripts")
        spe = .add_metrics(spe, "total_transcripts", "cell_level")
    }
    if ("total_genes" %in% metrics) {
        total_genes = Matrix::colSums(counts(spe) != 0)
        colData(spe)$total_genes = total_genes
        spe = .add_metrics(spe, "total_genes", "cell_level")
    }
    if ("total_cells" %in% metrics) {
        total_cells = Matrix::rowSums(counts(spe) != 0)
        SummarizedExperiment::rowData(spe)$total_cells = total_cells
        spe = .add_metrics(spe, "total_cells", "gene_level")
    }
    if ("meanExprsPct_cells" %in% metrics) {
        meanExprsPct_cells = Matrix::rowMeans(counts(spe) != 0)
        SummarizedExperiment::rowData(spe)$meanExprsPct_cells = meanExprsPct_cells
        spe = .add_metrics(spe, "meanExprsPct_cells", "gene_level")
    }
    return(spe)
}
spe = .initialise_CellSPA_list(spe)
spe = .add_dataset_metrics(spe)
spe = .cal_baseline(spe)

spe@metadata$CellSPA$method = tool_used

# Adding missing colums to spe colData
spe@colData['slide'] = 1
spe@colData['total_genes'] = colSums(spe@assays@data$counts != 0)
spe@colData['total_transcripts'] = colSums(spe@assays@data$counts)
spe@colData['total_reads'] = colSums(spe@assays@data$counts)
spe@colData['cell_id'] = colnames(spe@assays@data$counts)

# Assign segmentation output data only for cells present in metadata
coords = metadata_coord[c("x", "y", "cell_id")]
colnames(coords) = c("coord_x", "coord_y", "cell_id")
spe@metadata$CellSegOutput = coords

# Preprocess spe object
spe = processingSPE(spe, qc_range = list(total_transcripts = c(20, 2000), total_genes = c(20, Inf)))
#
# BASELINE METRICS COMPUTATION
#
spe = generatePolygon(spe, use_BPPARAM = BiocParallel::MulticoreParam(workers = parallel::detectCores() - 1))
spe = calBaselineAllMetrics(spe, verbose = TRUE)
#
# REFERENCE FULL ADD METADATA AND METRICS
ref_spe = SingleCellExperiment(assays = list(logcounts = logcounts(spe)))
#
ref_spe@metadata = data.frame(colData(spe))
#
ref_spe@metadata$celltype = factor(ref_spe@metadata$predicted.id)
#
ref_spe = processingRef(ref_spe, celltype = ref_spe@metadata$celltype , subset_row = rownames(spe))
#
metadata_dataframe = data.frame(celltype = colnames(ref_spe))
ref_spe@metadata = metadata_dataframe
#
# COMPUTE LARGE METRICS
# Calculating expression similarity
spe = calExpressionCorrelation(spe, ref_spe, ref_celltype = ref_spe@metadata$celltype, method = c("pearson", "cosine"), spe_exprs_values = "logcounts", ref_exprs_values = "mean")
spe = calExpressionCorrelation(spe, ref_spe, ref_celltype = ref_spe@metadata$celltype, method = c("pearson", "cosine"), spe_exprs_values = "logcounts", ref_exprs_values = "prop_detected")
spe = calAggExpressionCorrelation(spe, celltype = "mean_celltype_correlation", sce_ref = ref_spe, ref_celltype = "celltype", method = c("pearson"), spe_exprs_values = "logcounts", ref_exprs_values = "mean")
spe = calAggExpressionCorrelation(spe, celltype = "mean_celltype_correlation", sce_ref = ref_spe, ref_celltype = "celltype", method = c("pearson"), spe_exprs_values = "logcounts", ref_exprs_values = "prop_detected")
#
# Calculate marker F1 purity
# Positive
positive_marker_list = generateMarkerList(ref_spe, type = "positive")
spe = calMarkerPurity(spe, celltype = "mean_celltype_correlation", marker_list = positive_marker_list, marker_list_name = "positive")
# Negative
negative_marker_list = generateMarkerList(ref_spe, type = "negative", t = 1)
spe = calMarkerPurity(spe, celltype = "mean_celltype_correlation", marker_list = negative_marker_list, marker_list_name = "negative")
#
# Calculate marker expressed pct
# Positive
spe = calMarkerPct(spe, celltype = "mean_celltype_correlation", marker_list = positive_marker_list, marker_list_name = "positive")
# Negative
spe = calMarkerPct(spe, celltype = "mean_celltype_correlation", marker_list = negative_marker_list, marker_list_name = "negative")


# 6th
# Load required package and map cell_id from global image reference to each patch cell_id segmented on each tool (done in parallel)
process_paths = function(tool_used) {
    if (tool_used == "watershed"){
        #Breast
        output_path = path_to_watershed_segmentation_outcomes #Breast
        output_paths = list.files(output_path, full.names = TRUE)
        output_paths = output_paths[grep("spot2nucl_", output_paths)]
    }
    if (tool_used == "bidcell"){
        #Breast
        output_path = path_to_bidcell_segmentation_outcomes #Breast
        output_paths = list.files(output_path, full.names = TRUE)
    }
    if (tool_used == "cellpose"){
        #Breast
        output_path = path_to_cellpose_segmentation_outcomes #Breast
        output_paths = list.files(output_path, full.names = TRUE)
        output_paths = output_paths[grep("\\.txt$", output_paths)]
    }
    if (tool_used == "sam"){
        #Breast
        output_path = path_to_sam_segmentation_outcomes #Breast
        output_paths = list.files(output_path, full.names = TRUE)
    }
    if (tool_used == "sam2"){
        #Breast
        output_path = path_to_sam2_segmentation_outcomes #Breast
        output_paths = list.files(output_path, full.names = TRUE)
    }
    if (tool_used == "scs"){
        #Breast
        output_path = path_to_scs_segmentation_outcomes #Breast
        output_paths = list.files(output_path, full.names = TRUE)
        output_paths = output_paths[grep("spot2cell_", output_paths)]
    }
    if (tool_used == "scs_nucleus"){
        #Breast
        output_path = path_to_scs_nucleus_segmentation_outcomes #Breast
        output_paths = list.files(output_path, full.names = TRUE)
        output_paths = output_paths[grep("spot2nucl_", output_paths)]
    }
    if (tool_used == "baysor"){
        #Breast
        output_path = path_to_baysor_segmentation_outcomes #Breast
        output_paths = list.files(output_path, full.names = TRUE)
    }
    return(output_paths)
}

process_sample = function(n) {
    if (tool_used == "watershed"){
        sample_n = sub(".*/", "", output_paths[n])
        sample_n = sub(".*nucl_(.*?)\\..*", "\\1", sample_n)
        print(paste0('Working with watershed nucleus segmentation outputs'))
        print(paste0('Working in sample ', sample_n))
        output_paths_n = output_paths[n]
        #Load txt patch file
        output = read.csv(output_paths_n, header = FALSE, sep="\t")
        colnames(output) = c("pixel","cell_id")
        output = separate(output, "pixel", into = c("coord_x", "coord_y"), sep = ":")
        output = data.frame(lapply(output, as.double))         
    }
    if (tool_used == "bidcell"){
        sample_n = sub(".*/", "", output_paths[n])
        print(paste0('Working with bidcell segmentation outputs'))
        print(paste0('Working in sample ', sample_n))
        output_paths_n = list.files(output_paths[n], full.names = TRUE)[8]
        output_paths_n = list.files(output_paths_n, full.names = TRUE)
        output_paths_n = list.files(output_paths_n, full.names = TRUE)[3]
        output_paths_n_tiff = list.files(output_paths_n, full.names = TRUE)[3]
        #Load tiff patch file
        output = data.frame(readTiffOutput(output_paths_n_tiff))
    }
    if (tool_used == "cellpose"){
        sample_n = sub(".*/", "", output_paths[n])
        sample_n = sub(".*cell_(.*?)\\..*", "\\1", sample_n)
        print(paste0('Working with cellpose segmentation outputs'))
        print(paste0('Working in sample ', sample_n))
        output_paths_n = output_paths[n]
        #Load txt patch file
        output = read.csv(output_paths_n, header = FALSE, sep="\t")
        colnames(output) = c("pixel","cell_id")
        output = separate(output, "pixel", into = c("coord_x", "coord_y"), sep = ":")
        output = data.frame(lapply(output, as.double))
    }
    if (tool_used == "sam"){
        sample_n = sub(".*/", "", output_paths[n])
        print(paste0('Working with sam segmentation outputs'))
        print(paste0('Working in sample ', sample_n))
        output_paths_n = output_paths[n]
        #Load txt patch file
        output = read.csv(output_paths_n, header = FALSE, sep="\t")
        colnames(output) = c("pixel","cell_id")
        output = separate(output, "pixel", into = c("coord_x", "coord_y"), sep = ":")
        output = data.frame(lapply(output, as.double))
    }
    if (tool_used == "sam2"){
        sample_n = sub(".*/", "", output_paths[n])
        print(paste0('Working with sam2 segmentation outputs'))
        print(paste0('Working in sample ', sample_n))
        output_paths_n = output_paths[n]
        #Load txt patch file
        output = read.csv(output_paths_n, header = FALSE, sep="\t")
        colnames(output) = c("pixel","cell_id")
        output = separate(output, "pixel", into = c("coord_x", "coord_y"), sep = ":")
        output = data.frame(lapply(output, as.double))
    }
    if (tool_used == "scs"){
        sample_n = sub(".*/", "", output_paths[n])
        sample_n = sub(".*cell_(.*?)\\..*", "\\1", sample_n)
        print(paste0('Working with scs segmentation outputs'))
        print(paste0('Working in sample ', sample_n))
        output_paths_n = output_paths[n]
        #Load txt patch file
        output = read.csv(output_paths_n, header = FALSE, sep="\t")
        colnames(output) = c("pixel","cell_id")
        output = separate(output, "pixel", into = c("coord_x", "coord_y"), sep = ":")
        output = data.frame(lapply(output, as.double))       
    }
    if (tool_used == "scs_nucleus"){
        sample_n = sub(".*/", "", output_paths[n])
        sample_n = sub(".*nucl_(.*?)\\..*", "\\1", sample_n)
        print(paste0('Working with nucleus segmentation outputs'))
        print(paste0('Working in sample ', sample_n))
        output_paths_n = output_paths[n]
        #Load txt patch file
        output = read.csv(output_paths_n, header = FALSE, sep="\t")
        colnames(output) = c("pixel","cell_id")
        output = separate(output, "pixel", into = c("coord_x", "coord_y"), sep = ":")
        output = data.frame(lapply(output, as.double))       
    }
    if (tool_used == "baysor"){
        sample_n = sub(".*/", "", output_paths[n])
        print(paste0('Working with baysor segmentation outputs'))
        print(paste0('Working in sample ', sample_n))
        output_paths_n = output_paths[n]
        #Load txt patch file
        output = read.csv(output_paths_n, header = FALSE, sep="\t")
        colnames(output) = c("pixel","cell_id")
        output = separate(output, "pixel", into = c("coord_x", "coord_y"), sep = ":")
        output = data.frame(lapply(output, as.double))
    }
    #
    #Adapt output file
    y_ind = as.integer(sub(":.*", "", sample_n))
    if (tool_used == "scs" | tool_used == "scs_nucleus" | tool_used == "watershed"){
        x_ind = as.integer(sub(".*:(.*?):.*", "\\1", sample_n))
    } else {
        x_ind = as.integer(sub(".*:(.*?)_.*", "\\1", sample_n))
    }
    output["coord_x"] = output["coord_x"] + x_ind
    output["coord_y"] = output["coord_y"] + y_ind
    #
    #Adapt cell_names for each patch
    output$cell_id = paste0(output$cell_id, "_", n)
    #
    result = rbind(result, data.frame(output))
    #
    return(result)
}

# Bucle for each tool. Obtain new re-cuantified count matrix and the all-in-one segmentation outcomes for all pathes for each tool
for (tool in list_of_tools){
    mem_in_gb = mem_used() / (1024^3)
    print(mem_in_gb)
    # Set with tool is running
    tool_used = tool
    #
    # Get all patches for each tool
    output_paths = process_paths(tool_used = tool)
    #
    # Generate full segmentation outputs for each tool
    num_cores = detectCores() - 1
    print(num_cores)
    print('Start_mapping local cellID with grobal cellID')
    # Initialize object for store each patch for each single tool
    result = data.frame()
    #
    combined_output = do.call(rbind, mclapply(1:length(output_paths), process_sample, mc.cores = num_cores))
    print(length(unique(combined_output[,3])))
    #
    library(dplyr)
    library(tidyr)
    #
    # Step 1: Merge the spatial coordinates and combined output dataframes
    spatial_coords=df[c("x","y","geneID","MIDCounts")]
    colnames(spatial_coords) = c("coord_x","coord_y","geneID","MIDCounts")
    combined_output$coord_x = as.integer(combined_output$coord_x)
    combined_output$coord_y= as.integer(combined_output$coord_y)
    merged_data = inner_join(spatial_coords, combined_output, by = c("coord_x", "coord_y"))
    print (dim(merged_data))
    # Step 2: Group by 'cellid' and 'gene', then sum the 'count' values
    cell_gene_counts = merged_data %>%
                        group_by(cell_id, geneID) %>%
                        summarise(counts = sum(MIDCounts)) %>%
                        ungroup()

    # Step 3: Pivot the data to create a gene expression matrix with 'cellid' as rows and 'gene' as columns
    new_count_matrix = cell_gene_counts %>%
                       pivot_wider(names_from = geneID, values_from = counts, values_fill = list(counts = 0))
    new_count_matrix = data.frame(new_count_matrix) #this can be mem optimized by removing transforming into data.frame
    rownames(new_count_matrix) = new_count_matrix[,1]
    new_count_matrix= new_count_matrix[,-1]
    new_count_matrix = t(new_count_matrix) 
    print(dim(new_count_matrix))
    #
    # Subset segmentation outcomes base only on those cells that matched and contains transcriptomic information
    combined_output = combined_output[combined_output$cell_id %in% colnames(new_count_matrix),]
    #
    # Generate objects for each tool files. Are, new_count_matrix(re-quantified segmentation count matrix) and combined_outcome(segmentation)
    current_tool_new_count_matrix = paste0("count_matrix_", tool)
    assign_expr =  paste0(current_tool_new_count_matrix, " = new_count_matrix")
    eval(parse(text = assign_expr))
    current_tool_new_name = paste0("combined_output_", tool)
    assign_expr =  paste0(current_tool_new_name, " = combined_output")
    eval(parse(text = assign_expr))
}

#Save files into path and create directory
wd = path_output_files
dataset_folder_name = tissue_subpath

output_path = paste0(wd, dataset_folder_name)

#Create root path
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# SAVEPOINT
save.image(paste0(output_path, "/savepoint_1.RData"))
#load(paste0(output_path, "/savepoint_1.RData"))

###
###
#List of tools
list_of_tools = list_of_segmentation_tools

# Loop through the list_of_tools desired
for (tool in list_of_tools) {
    #Check mem used
    mem_in_gb = mem_used() / (1024^3)
    print(mem_in_gb)
    #
    tool_used = tool
    # START WITH CELLSPA 
    print(paste0("WORKING IN COMBINED PATCHES TOOL: ", paste0("combined_output_", tool_used)))
    segmentation_output = get(paste0("combined_output_", tool_used)) #tif_output
    #
    # Remove duplicates from new count matrix and generate a new spatialExperiment object with it
    obj_name = paste0(tool_used, "_spe")
    matrix_obj_name = paste0("count_matrix_", tool_used)
    assign(matrix_obj_name, get(matrix_obj_name)[, !is.na(colnames(get(matrix_obj_name)))])
    assign(obj_name, SpatialExperiment::SpatialExperiment(assay = list(counts = get(matrix_obj_name)))) 
    #
    # Compute logcounts for object
    assign_expr = paste0(obj_name, "@assays@data$logcounts = log2(get(obj_name)@assays@data$counts + 1)")
    eval(parse(text = assign_expr))
    #
    # Assign ssegmentation output data only for cells present in metadata
    assign_expr = paste0(obj_name, "@metadata$CellSegOutput = segmentation_output")
    eval(parse(text = assign_expr))
    #
    # Add extra data required to object
    assign(obj_name, .initialise_CellSPA_list(get(obj_name))) 
    assign(obj_name, .add_dataset_metrics(get(obj_name))) 
    assign(obj_name, .cal_baseline(get(obj_name))) 
    #
    assign_expr = paste0(obj_name, "@metadata$CellSPA$method = tool_used")
    eval(parse(text = assign_expr))
    #
    # Adding missing colums to spe colData
    assign_expr = paste0(obj_name, "@colData$slide = 1")
    eval(parse(text = assign_expr))
    assign_expr = paste0(obj_name, "@colData$total_genes = colSums(get(obj_name)@assays@data$counts != 0)")
    eval(parse(text = assign_expr))
    assign_expr = paste0(obj_name, "@colData$total_transcripts = colSums(get(obj_name)@assays@data$counts)")
    eval(parse(text = assign_expr))
    assign_expr = paste0(obj_name, "@colData$total_reads = colSums(get(obj_name)@assays@data$counts)")
    eval(parse(text = assign_expr))
    assign_expr = paste0(obj_name, "@colData$cell_id = colnames(get(obj_name)@assays@data$counts)")
    eval(parse(text = assign_expr))
    #
    # Process the data, filtering low quality cells
    assign(obj_name, processingSPE(get(obj_name), qc_range = list(total_transcripts = c(20, 2000), total_genes = c(20, Inf))))
    #
    # BASELINE METRICS COMPUTATION
    # Fix any possible error before compute the metrics
    assign_expr = paste0(obj_name, "@metadata$CellSegOutput = na.omit(get(obj_name)@metadata$CellSegOutput)")
    eval(parse(text = assign_expr))
    assign_expr = paste0(obj_name, "= get(obj_name)[, unique(get(obj_name)@metadata$CellSegOutput$cell_id)]")
    eval(parse(text = assign_expr))
    assign_expr = paste0(obj_name, "@metadata$CellSegOutput= unique(get(obj_name)@metadata$CellSegOutput)")
    eval(parse(text = assign_expr))
    #
    # Filter unrealistic extreme big segmentations (makes polygon computation never ends)
    print(dim(get(obj_name)@metadata$CellSegOutput))
    unreal_segmentations = names(table(get(obj_name)@metadata$CellSegOutput$cell_id)[table(get(obj_name)@metadata$CellSegOutput$cell_id) > 5000])
    assign_expr = paste0(obj_name, "@metadata$CellSegOutput = get(obj_name)@metadata$CellSegOutput[!(get(obj_name)@metadata$CellSegOutput$cell_id %in% unreal_segmentations), ]")
    eval(parse(text = assign_expr))
    print(dim(get(obj_name)@metadata$CellSegOutput))
    #
    assign(obj_name, generatePolygon(get(obj_name), use_BPPARAM = BiocParallel::MulticoreParam(workers = 5))) #parallel::detectCores() - 1
    assign(obj_name, calBaselineAllMetrics(get(obj_name), verbose = TRUE))
    #
    # REFERENCE FULL ADD METADATA AND METRICS
    ref_obj_name = paste0(tool, "_ref_data")
    assign(ref_obj_name, SingleCellExperiment(assays = list(logcounts = logcounts(spe)))) #spe is the full, same filtered reference data
    #
    assign_expr = paste0(ref_obj_name, "@metadata = data.frame(colData(spe))")
    eval(parse(text = assign_expr))
    #
    assign_expr = paste0(ref_obj_name, "@metadata$celltype = factor(get(ref_obj_name)@metadata$predicted.id)") #annotations
    eval(parse(text = assign_expr))
    #
    assign(ref_obj_name, processingRef(get(ref_obj_name), celltype = get(ref_obj_name)@metadata$celltype , subset_row = rownames(spe))) 
    #
    metadata_dataframe = data.frame(celltype = colnames(get(ref_obj_name)))
    assign_expr = paste0(ref_obj_name, "@metadata = metadata_dataframe")
    eval(parse(text = assign_expr))
    #
    # COMPUTE LARGE METRICS
    # Calculating expression similarity
    assign(obj_name, calExpressionCorrelation(get(obj_name), get(ref_obj_name), ref_celltype = get(ref_obj_name)@metadata$celltype, method = c("pearson", "cosine"), spe_exprs_values = "logcounts", ref_exprs_values = "mean"))
    assign(obj_name, calExpressionCorrelation(get(obj_name), get(ref_obj_name), ref_celltype = get(ref_obj_name)@metadata$celltype, method = c("pearson", "cosine"), spe_exprs_values = "logcounts", ref_exprs_values = "prop_detected"))
    assign(obj_name, calAggExpressionCorrelation(get(obj_name), celltype = "mean_celltype_correlation", sce_ref = get(ref_obj_name), ref_celltype = "celltype", method = c("pearson"), spe_exprs_values = "logcounts", ref_exprs_values = "mean"))
    assign(obj_name, calAggExpressionCorrelation(get(obj_name), celltype = "mean_celltype_correlation", sce_ref = get(ref_obj_name), ref_celltype = "celltype", method = c("pearson"), spe_exprs_values = "logcounts", ref_exprs_values = "prop_detected"))
    #
    # Calculate marker F1 purity
    # Positive
    positive_marker_list = generateMarkerList(get(ref_obj_name), type = "positive")
    assign(obj_name, calMarkerPurity(get(obj_name), celltype = "mean_celltype_correlation", marker_list = positive_marker_list, marker_list_name = "positive"))
    # Negative
    negative_marker_list = generateMarkerList(get(ref_obj_name), type = "negative", t = 1)
    assign(obj_name, calMarkerPurity(get(obj_name), celltype = "mean_celltype_correlation", marker_list = negative_marker_list, marker_list_name = "negative"))
    #
    # Calculate marker expressed pct
    # Positive
    assign(obj_name, calMarkerPct(get(obj_name), celltype = "mean_celltype_correlation", marker_list = positive_marker_list, marker_list_name = "positive"))
    #Negative
    assign(obj_name, calMarkerPct(get(obj_name), celltype = "mean_celltype_correlation", marker_list = negative_marker_list, marker_list_name = "negative"))
}

# SAVEPOINT
save.image(paste0(output_path, "/savepoint_2.RData"))


###
###
## WORKING ON VISUALIZATION PART
#load(paste0(output_path, "/savepoint_2.RData"))

list_of_tools = list_of_segmentation_tools

# PLOTS PART 1 -OVERALL plots and Gene-level QC metrics-

############
#FIRST PAGE
# 1st plot
# Plot the number of cell segmented by each method that are mapped to the reference as cells 
# Prepare data for plot
obj_name_list_of_tools_ncells = paste0("list_of_tools_ncells")
list_of_tools_ncells = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")
    assign_expr =  paste0(obj_name_list_of_tools_ncells, " = append(list_of_tools_ncells, nrow(colData(get(obj_name))))")
    eval(parse(text = assign_expr))
}

ncell_comp = data.frame(tool = list_of_tools, 
                        cells = list_of_tools_ncells) 

# Ordering 
ncell_comp$tool = factor(ncell_comp$tool, levels = ncell_comp$tool)

# Assigning default and different colors to bar plot 
n_cells_plot = ggplot(data=ncell_comp, aes(x=tool, y=cells,fill=tool))+ 
                    geom_bar(stat="identity") +
                    labs(title = "Number of cells", x = "", y = "") +
                    theme(plot.title = element_text(hjust = 0.5))

# Create path for figures for corresponding dataset
# Create tools paths
if (!dir.exists(paste0(output_path, "/images"))) {
    dir.create(paste0(output_path, "/images"), recursive = TRUE)
    cat("Created directory:", path, "\n")
}

# Save plot
pdf(paste0(output_path, "/images/Ncells_tools_comparison.pdf"))
n_cells_plot 
dev.off()

# 2nd plot
# Plot the % of transcripts assigned by each method related to the global transcipts identified
# Prepare data for plot
# Reference total transcripts assigned over all reference cells
referece_total_transcripts = sum(colData(spe)$total_transcripts)

obj_name_tool_total_transcripts = paste0("list_tool_total_transcripts")
list_tool_total_transcripts = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")
    tool_total_transcripts = sum(colData(get(obj_name))$total_transcripts)
    assign_expr =  paste0(obj_name_tool_total_transcripts, " = append(list_tool_total_transcripts, tool_total_transcripts)")
    eval(parse(text = assign_expr))
}

# Compute the ratio for each tool respect to reference
list_tool_total_transcripts_pct = (list_tool_total_transcripts/referece_total_transcripts)*100

transcripts_pct_comp = data.frame(tool = list_of_tools, 
                        pct = list_tool_total_transcripts_pct) 

# Ordering 
transcripts_pct_comp$tool = factor(transcripts_pct_comp$tool, levels = transcripts_pct_comp$tool)

# Assigning default and different colors to bar plot 
transcripts_assigned_plot = ggplot(data=transcripts_pct_comp, aes(x=tool, y=pct,fill=tool))+ 
                                geom_bar(stat="identity") +
                                labs(title = "% of transcripts assigned", x = "", y = "") +
                                theme(plot.title = element_text(hjust = 0.5))

# Save plot
pdf(paste0(output_path, "/images/Tran_assig_pct_tools_comparison.pdf"))
transcripts_assigned_plot 
dev.off()

# 3rd plot
# Plot the % of transcripts assigned by each method related to the global transcipts identified
# Prepare data for plot
obj_name_tool_total_transcripts = paste0("list_tool_total_transcripts")
list_tool_total_transcripts = c()
obj_name_tool_name_total_transcripts = paste0("list_tool_name_total_transcripts")
list_tool_name_total_transcripts = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")
    tool_total_transcripts = as.numeric(colData(get(obj_name))$total_transcripts)
    assign_expr =  paste0(obj_name_tool_total_transcripts, " = append(list_tool_total_transcripts, tool_total_transcripts)")
    eval(parse(text = assign_expr))
    # Create tool name vector for correspoding set of cells
    tool_name_total_transcripts = rep(tool, length(tool_total_transcripts))
    assign_expr =  paste0(obj_name_tool_name_total_transcripts, " = append(list_tool_name_total_transcripts, tool_name_total_transcripts)")
    eval(parse(text = assign_expr))
}

transcripts_total_comp = data.frame(tool = list_tool_name_total_transcripts, 
                        total = list_tool_total_transcripts)

# Ordering 
transcripts_total_comp$tool = factor(transcripts_total_comp$tool , levels = list_of_tools)

# Assigning default and different colors to bar plot 
total_transcripts_plot = ggplot(data=transcripts_total_comp, aes(x=tool, y=total, fill=tool))+ 
                    geom_boxplot() +
                    labs(title = "Number of total transcripts per cell", x = "", y = "") +
                    theme(plot.title = element_text(hjust = 0.5))

# Save plot
pdf(paste0(output_path, "/images/Tran_assig_total_tools_comparison.pdf"))
total_transcripts_plot 
dev.off()

# 4th plot
# Plot the % of transcripts assigned by each method related to the global transcipts identified
# Prepare data for plot
obj_name_tool_total_genes = paste0("list_tool_total_genes")
list_tool_total_genes = c()
obj_name_tool_name_total_genes = paste0("list_tool_name_total_genes")
list_tool_name_total_genes = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")
    tool_total_genes = as.numeric(colData(get(obj_name))$total_genes)
    assign_expr =  paste0(obj_name_tool_total_genes, " = append(list_tool_total_genes, tool_total_genes)")
    eval(parse(text = assign_expr))
    # Create tool name vector for correspoding set of cells
    tool_name_total_genes = rep(tool, length(tool_total_genes))
    assign_expr =  paste0(obj_name_tool_name_total_genes, " = append(list_tool_name_total_genes, tool_name_total_genes)")
    eval(parse(text = assign_expr))
}

genes_total_comp = data.frame(tool = list_tool_name_total_genes, 
                        total = list_tool_total_genes)

# Ordering 
genes_total_comp$tool = factor(genes_total_comp$tool , levels = list_of_tools)

# Assigning default and different colors to bar plot 
total_genes_plot = ggplot(data=genes_total_comp, aes(x=tool, y=total, fill=tool))+ 
                    geom_boxplot() +
                    labs(title = "Number of total genes per cell", x = "", y = "") +
                    theme(plot.title = element_text(hjust = 0.5))

# Save plot
pdf(paste0(output_path, "/images/Genes_assig_total_tools_comparison.pdf"))
total_genes_plot 
dev.off()

# 5th plot
# Plot gene-level qc metrics comparing segmentation and nucleus % of cells expressing each gene
# Prepare data for plot

# Nucleus reference
ref_nucleus = "watershed"

obj_name_tool_total_meanExprsPct_cells = paste0("list_tool_total_meanExprsPct_cells")
list_tool_total_meanExprsPct_cells = c()
obj_name_tool_name_total_meanExprsPct_cells = paste0("list_tool_name_total_meanExprsPct_cells")
list_tool_name_total_meanExprsPct_cells = c()
obj_name_referece_genes_assig_nucleus = paste0("list_referece_genes_assig_nucleus")
list_referece_genes_assig_nucleus = c()

# Remove reference from list of tools to be used here
list_of_tools_minus_ref = list_of_tools[!list_of_tools %in% ref_nucleus]

for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")
    tool_meanExprsPct_cells = (as.numeric(rowData(get(obj_name))$meanExprsPct_cells))*100
    assign_expr =  paste0(obj_name_tool_total_meanExprsPct_cells, " = append(list_tool_total_meanExprsPct_cells, tool_meanExprsPct_cells)")
    eval(parse(text = assign_expr))

    obj_name_ref_nucleus = paste0(ref_nucleus, "_spe") 
    tool_meanExprsPct_cells_gene_names = names(rowData(get(obj_name))$meanExprsPct_cells)
    referece_genes_assig_nucleus = (as.numeric(rowData(get(obj_name_ref_nucleus))$meanExprsPct_cells[tool_meanExprsPct_cells_gene_names]))*100
    assign_expr =  paste0(obj_name_referece_genes_assig_nucleus, " = append(list_referece_genes_assig_nucleus, referece_genes_assig_nucleus)")
    eval(parse(text = assign_expr))

    # Create tool name vector for correspoding set of cells
    tool_name_meanExprsPct_cells = rep(tool, length(tool_meanExprsPct_cells))
    assign_expr =  paste0(obj_name_tool_name_total_meanExprsPct_cells, " = append(list_tool_name_total_meanExprsPct_cells, tool_name_meanExprsPct_cells)")
    eval(parse(text = assign_expr))

    # ADD gene names only values for reference

}

pct_assigned_genes = data.frame(tool = list_tool_name_total_meanExprsPct_cells, 
                                ref_pct_assig_genes = list_referece_genes_assig_nucleus,
                                pct_assig_genes = list_tool_total_meanExprsPct_cells)

# Ordering 
pct_assigned_genes$tool = factor(pct_assigned_genes$tool , levels = list_of_tools)

# Assigning default and different colors to bar plot 
pct_assig_genes = ggplot(data=pct_assigned_genes, aes(x=ref_pct_assig_genes, y=pct_assig_genes, color=tool)) +
                    geom_point() +
                    geom_abline(intercept = 0, slope = 1, color = "red") +
                    facet_wrap(~tool, ncol = length(list_of_tools)) +
                    labs(title = "% of cells expressing each gene", x = "Nuclei", y = "Segmented cells") +
                    coord_fixed(ratio = 1)  

# Save plot
pdf(paste0(output_path, "/images/pct_Assigned_genes_tools_comparison.pdf"))
pct_assig_genes 
dev.off()

# 6th plot and 7th plot
# Plot cell morphology metrics - elongation between nucleus and segmented cells by cell types
# Prepare data for plot
# Nucleus reference
ref_nucleus = "watershed"
obj_name_ref_nucleus = paste0(ref_nucleus, "_spe") 
referece_elongation = data.frame(colData(get(obj_name_ref_nucleus))) %>%
                                            group_by(mean_celltype_correlation) %>%
                                            summarize(average_value = mean(elongation, na.rm = TRUE))
referece_elongation_names = referece_elongation$mean_celltype_correlation
referece_elongation = referece_elongation$average_value
names(referece_elongation) = referece_elongation_names

obj_name_tool_elongation_cells = paste0("list_tool_total_elongation_cells")
list_tool_total_elongation_cells = c()
obj_name_tool_celltype_cells = paste0("list_tool_celltype_cells")
list_tool_celltype_cells = c()
obj_name_tool_name_elongation_cells = paste0("list_tool_name_elongation_cells")
list_tool_name_elongation_cells = c()
obj_name_tool_total_transcripts = paste0("list_tool_total_transcripts")
list_tool_total_transcripts = c()
obj_reference_elongation = paste0("list_reference_elongation")
list_reference_elongation = c()

# Remove reference from list of tools to be used here
list_of_tools_minus_ref = list_of_tools[!list_of_tools %in% ref_nucleus]

# Initialize empty correlation
correlation_vec = c()
tools_average_total_transcripts = c()

for (tool in list_of_tools_minus_ref){
    print(tool)
    obj_name  = paste0(tool, "_spe")

    tool_meanExprsPct_cells_matrix = data.frame(colData(get(obj_name))) %>%
                                            group_by(mean_celltype_correlation) %>%
                                            summarize(average_value = mean(elongation, na.rm = TRUE))
    tool_meanExprsPct_cells_names = tool_meanExprsPct_cells_matrix$mean_celltype_correlation
    tool_meanExprsPct_cells = tool_meanExprsPct_cells_matrix$average_value
    names(tool_meanExprsPct_cells) = tool_meanExprsPct_cells_names

    # Select common celltypes
    overlap_names = tool_meanExprsPct_cells_names[tool_meanExprsPct_cells_names %in% referece_elongation_names]
    tool_meanExprsPct_cells = tool_meanExprsPct_cells[overlap_names]
    referece_elongation_iter = referece_elongation[overlap_names]

    tool_celltype_cor_value = cor(referece_elongation_iter, tool_meanExprsPct_cells, method="spearman")

    correlation_vec = append(correlation_vec, tool_celltype_cor_value)

    # Create tool name vector for correspoding set of cells
    tool_name_meanExprsPct_cells = rep(tool, length(tool_meanExprsPct_cells))
    assign_expr =  paste0(obj_name_tool_name_elongation_cells, " = append(list_tool_name_elongation_cells, tool_name_meanExprsPct_cells)")
    eval(parse(text = assign_expr))
    assign_expr =  paste0(obj_name_tool_elongation_cells, " = append(list_tool_total_elongation_cells, tool_meanExprsPct_cells)")
    eval(parse(text = assign_expr))
    assign_expr =  paste0(obj_name_tool_celltype_cells, " = append(list_tool_celltype_cells, overlap_names)")
    eval(parse(text = assign_expr))
    assign_expr =  paste0(obj_reference_elongation, " = append(list_reference_elongation, referece_elongation_iter)")
    eval(parse(text = assign_expr))

    # Compute mean total transcript for each tool
    tool_total_transcripts = mean(as.numeric(colData(get(obj_name))$total_transcripts))
    assign_expr =  paste0(obj_name_tool_total_transcripts, " = append(list_tool_total_transcripts, tool_total_transcripts)")
    eval(parse(text = assign_expr))
    tools_average_total_transcripts = append(tools_average_total_transcripts, tool_total_transcripts)
}

# Matrix for elongation morphology
elongations = data.frame(tool = list_tool_name_elongation_cells,
                                celltype = list_tool_celltype_cells,
                                ref_elongation = list_reference_elongation,
                                segmentation_elongations = list_tool_total_elongation_cells)

names(correlation_vec) = list_of_tools_minus_ref
elongations = elongations %>%
  mutate(tool_correlation = recode(tool, !!!correlation_vec))

# Ordering 
elongations$tool = factor(elongations$tool , levels = list_of_tools)

# Matrix for Average total transcript comparison with celltype elongation morphology correlarion
transcipts_total_comp = data.frame(tool_correlations = correlation_vec, 
                                   tool_average_transcripts = tools_average_total_transcripts)

# Ordering 
transcipts_total_comp$tool = factor(rownames(transcipts_total_comp) , levels = list_of_tools)

# Assigning default and different colors to bar plot 
segmentation_enlongations = ggplot(data=elongations, aes(x=ref_elongation, y=segmentation_elongations, color=celltype)) +
                    geom_point() +
                    facet_wrap(~tool, ncol = length(list_of_tools_minus_ref)) +
                    labs(title = "Elongation", x = "Nuclei", y = "Segmented cells") +
                    coord_fixed(ratio = 1) +
                    geom_text(data = elongations, aes(x = Inf, y = Inf, 
                                label = paste("Corr:", round(tool_correlation, 2))), 
                                hjust = 1.1, vjust = 1.5, size = 2.5, color = "black", inherit.aes = FALSE) +
                    theme(legend.position = "right") +
                    guides(color = guide_legend(ncol = 1)) +
                    theme(
                        legend.key.size = unit(0.5, "lines"),  # Smaller legend keys
                        legend.text = element_text(size = 8),
                        axis.text.x = element_text(size = 7),
                        axis.text.y = element_text(size = 7)   # Smaller legend text
                    )

correlation_elongations = ggplot(data=transcipts_total_comp, aes(x=tool_average_transcripts, y=tool_correlations, color=tool))+ 
                    geom_point(size = 8) +
                    labs(title = "", x = "Average of total trnascripts per cell", y = "Pearson correlation of cell-type elongations") +
                    theme(plot.title = element_text(hjust = 0.5)) +
                    geom_text(aes(label = tool), nudge_x = 0.1, nudge_y = 0.01) + 
                    theme(legend.position = "none")

# Save plot
pdf(paste0(output_path, "/images/elongatiopn_morphology_tools_comparison.pdf"))
segmentation_enlongations 
dev.off()

# Save plot
pdf(paste0(output_path, "/images/elongation_corr_vs_average_transcripts_tools_comparison.pdf"))
correlation_elongations 
dev.off()

######################
# Grid for PAGE 1 plots
# Plots to generate PAGE 1
p1 = n_cells_plot
p2 = transcripts_assigned_plot
p3 = total_transcripts_plot
p4 = total_genes_plot
p5 = pct_assig_genes
p6 = segmentation_enlongations
p7 = correlation_elongations

p1 = p1 + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) 
p2 = p2 + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
p3 = p3 + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
p4 = p4 + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
p5 = p5 + theme(legend.position = "none")

p6 = p6
p7 = p7 

# Arrange plots in a grid
pdf(paste0(output_path, "/images/PAGE_1.pdf"), width = 15, height = 15)
grid.arrange(
             arrangeGrob(p1, p2, p3, p4,  ncol = 4, nrow = 1),
             arrangeGrob(p5, ncol = 1, nrow = 1),
             arrangeGrob(p6, p7, ncol = 2, nrow = 1),
             nrow = 4
)
dev.off()


############
############
# SECOND PAGE
# EXPRESSION AND CELL COMPOSITION SIMILARITY WITH SCRNA-SEQ
# 1st plot
# Heatmap for average expression comparison between segmented and original nucleous data
# Prepare data for plot

avg_diag_value_all = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")
    #
    # Aggregated average cell type expression correlation by cell_types
    agg_mean_correlation = get(obj_name)@metadata$CellSPA$similarity_metrics$agg_mean_correlation
    agg_mean_correlation = agg_mean_correlation[, colSums(agg_mean_correlation != 0) > 0]
    agg_mean_correlation = agg_mean_correlation[colnames(agg_mean_correlation),]
    #
    # ADD sizes to rownames
    fontsizes = rep(10, nrow(agg_mean_correlation))
    #
    # Create text annotation object for displaying row names
    rowAnno = rowAnnotation(rows = anno_text(rownames(agg_mean_correlation), gp = gpar(fontsize = fontsizes)))
    columnAnno = columnAnnotation(rows = anno_text(colnames(agg_mean_correlation), gp = gpar(fontsize = fontsizes)))
    # Create color base palette 
    col_fun = colorRamp2(c(min(agg_mean_correlation), max(agg_mean_correlation)), c("white", "purple"))
    #
    # Diagonal average value vector for each tool
    avg_diag_value_all = append(avg_diag_value_all, mean(diag(agg_mean_correlation)))
    #
    # Pearson correlation between cell types
    celltype_corr_heatmap = Heatmap(agg_mean_correlation,
                                    name                         = "Pearson correlation",
                                    col                          = col_fun,
                                    show_row_names               = FALSE,
                                    show_column_names            = TRUE,
                                    row_names_gp                 = gpar(fontsize = 6),
                                    row_title_rot                = 0,
                                    cluster_rows                 = FALSE,
                                    cluster_row_slices           = FALSE,
                                    cluster_columns              = FALSE,
                                    row_dend_reorder             = FALSE,
                                    right_annotation             = rowAnno,
                                    row_order                    = rownames(agg_mean_correlation))
    #
    pdf(paste0(output_path, "/images/", "Average_corr_exp_", tool,"_heatmap.pdf"), width = 15, height = 15)
    print(celltype_corr_heatmap)
    dev.off()
}

############
# 2nd
# Plot average trsanscripts per cell agains pearson corr of celltypes correlation agains reference celltype data expression values
obj_tool_corr_cells = paste0("list_tool_total_corr_cells")
list_tool_total_corr_cells = c()
obj_name_tool_coor_cells = paste0("list_tool_names_celltypes_cells")
list_tool_names_celltypes_cells = c()
obj_name_tool_name_corr_cells = paste0("list_tool_names_corr_cells")
list_tool_names_corr_cells = c()
obj_tool_total_transcripts = paste0("list_tool_total_transcripts")
list_tool_total_transcripts = c()

tools_average_total_transcripts = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")

    tool_celltype_corr_cells_matrix = data.frame(colData(get(obj_name))) %>%
                                            group_by(mean_celltype_correlation) %>%
                                            summarize(average_value = mean(mean_cor_correlation, na.rm = TRUE))
    tool_corr_cells = tool_celltype_corr_cells_matrix$average_value
    tool_celltype_corr_cells = tool_celltype_corr_cells_matrix$mean_celltype_correlation
    assign_expr =  paste0(obj_tool_corr_cells, " = append(list_tool_total_corr_cells, tool_corr_cells)")
    eval(parse(text = assign_expr))
    assign_expr =  paste0(obj_name_tool_coor_cells, " = append(list_tool_names_celltypes_cells, tool_celltype_corr_cells)")
    eval(parse(text = assign_expr))

    # Create tool name vector for correspoding set of cells
    tool_name_meanExprsPct_cells = rep(tool, length(tool_corr_cells))
    assign_expr =  paste0(obj_name_tool_name_corr_cells, " = append(list_tool_names_corr_cells, tool_name_meanExprsPct_cells)")
    eval(parse(text = assign_expr))


    # Compute mean total transcript for each tool
    tool_total_transcripts = mean(as.numeric(colData(get(obj_name))$total_transcripts))
    assign_expr =  paste0(obj_tool_total_transcripts, " = append(list_tool_total_transcripts, tool_total_transcripts)")
    eval(parse(text = assign_expr))
    tools_average_total_transcripts = append(tools_average_total_transcripts, tool_total_transcripts) 
}

# Matrix for celtype corr
celltype_correlations = data.frame(tool = list_tool_names_corr_cells,
                                   celltype = list_tool_names_celltypes_cells,
                                   correlations = list_tool_total_corr_cells)

# Ordering 
celltype_correlations$tool = factor(celltype_correlations$tool , levels = list_of_tools)

celltype_correlations = celltype_correlations %>%
  group_by(tool) %>%
  summarize(sum_value = mean(correlations))

celltype_correlations$tools_average_total_transcripts = tools_average_total_transcripts
celltype_correlations = data.frame(celltype_correlations)

correlation_expression = ggplot(data=celltype_correlations, aes(x=tools_average_total_transcripts, y=sum_value, color=tool))+ 
                    geom_point(size = 8) +
                    labs(title = "", x = "Average of total trnascripts per cell", y = "Pearson correlation with reference") +
                    theme(plot.title = element_text(hjust = 0.5)) +
                    geom_text(aes(label = tool), nudge_x = 8, nudge_y = 0.01) + 
                    theme(legend.position = "none")

pdf(paste0(output_path, "/images/", "Expression_celltype_corr_vs_ref.pdf"))
print(correlation_expression)
dev.off()

############
# 3rd 
# Plot average trsanscripts per cell agains pearson corr of celltypes correlation agains reference celltype data proportion values
obj_tool_corr_cells = paste0("list_tool_total_corr_cells")
list_tool_total_corr_cells = c()
obj_name_tool_coor_cells = paste0("list_tool_names_celltypes_cells")
list_tool_names_celltypes_cells = c()
obj_name_tool_name_corr_cells = paste0("list_tool_names_corr_cells")
list_tool_names_corr_cells = c()
obj_tool_total_transcripts = paste0("list_tool_total_transcripts")
list_tool_total_transcripts = c()

tools_average_total_transcripts = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")

    tool_celltype_corr_cells_matrix = data.frame(colData(get(obj_name))) %>%
                                            group_by(mean_celltype_correlation) %>%
                                            summarize(average_value = mean(prop_detected_cor_correlation, na.rm = TRUE)) 
    tool_corr_cells = tool_celltype_corr_cells_matrix$average_value
    tool_celltype_corr_cells = tool_celltype_corr_cells_matrix$mean_celltype_correlation
    assign_expr =  paste0(obj_tool_corr_cells, " = append(list_tool_total_corr_cells, tool_corr_cells)")
    eval(parse(text = assign_expr))
    assign_expr =  paste0(obj_name_tool_coor_cells, " = append(list_tool_names_celltypes_cells, tool_celltype_corr_cells)")
    eval(parse(text = assign_expr))

    # Create tool name vector for correspoding set of cells
    tool_name_meanExprsPct_cells = rep(tool, length(tool_corr_cells))
    assign_expr =  paste0(obj_name_tool_name_corr_cells, " = append(list_tool_names_corr_cells, tool_name_meanExprsPct_cells)")
    eval(parse(text = assign_expr))


    # Compute mean total transcript for each tool
    tool_total_transcripts = mean(as.numeric(colData(get(obj_name))$total_transcripts))
    assign_expr =  paste0(obj_tool_total_transcripts, " = append(list_tool_total_transcripts, tool_total_transcripts)")
    eval(parse(text = assign_expr))
    tools_average_total_transcripts = append(tools_average_total_transcripts, tool_total_transcripts) 
}

# Matrix for celtype corr
celltype_correlations2 = data.frame(tool = list_tool_names_corr_cells,
                                   celltype = list_tool_names_celltypes_cells,
                                   correlations = list_tool_total_corr_cells)

# Ordering 
celltype_correlations2$tool = factor(celltype_correlations2$tool , levels = list_of_tools)

celltype_correlations2 = celltype_correlations2 %>%
  group_by(tool) %>%
  summarize(sum_value = mean(correlations))

celltype_correlations2$tools_average_total_transcripts = tools_average_total_transcripts
celltype_correlations2 = data.frame(celltype_correlations2)

correlation_proportions = ggplot(data=celltype_correlations2, aes(x=tools_average_total_transcripts, y=sum_value, color=tool))+ 
                    geom_point(size = 8) +
                    labs(title = "", x = "Average of total trnascripts per cell", y = "Pearson correlation with reference") +
                    theme(plot.title = element_text(hjust = 0.5)) +
                    geom_text(aes(label = tool), nudge_x = 8, nudge_y = 0.01) + 
                    theme(legend.position = "none")

pdf(paste0(output_path, "/images/", "Proportions_celltype_corr_vs_ref.pdf"))
print(correlation_proportions)
dev.off()

############
# 4th 
# Celltype proportion comparison
# Celltype frequency tool
tool_name_rep_total = c()
vector1_common_all = c()
vector2_common_all = c()
common_names_all = c()
tools_coor = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")
    ref_obj_name = paste0(tool, "_ref_data")

    tool_celltype_prop = (table(colData(get(obj_name))$mean_celltype_correlation)/dim(get(obj_name))[2])*100

    # Celltype frequency reference
    ref_celltype_names = colData(get(ref_obj_name))$celltype
    ref_celltype_prop = colData(get(ref_obj_name))$Freq *100
    names(ref_celltype_prop) = ref_celltype_names

    # Find the overlapping names
    common_names = intersect(names(ref_celltype_prop), names(tool_celltype_prop))
    common_names = factor(common_names, levels = common_names)

    # Subset the vectors based on the overlapping names
    vector1_common = ref_celltype_prop[common_names]
    vector2_common = tool_celltype_prop[common_names]

    vector1_common_all = append(vector1_common_all, as.numeric(vector1_common))
    vector2_common_all = append(vector2_common_all, as.numeric(vector2_common))

    common_names_all = append(common_names_all, common_names)

    # Cor value to add to plot
    tool_coor = cor(vector1_common,vector2_common)
    tool_coor = rep(tool_coor, length(common_names))
    tools_coor = append(tools_coor, tool_coor)

    # Tool name rep 
    tool_name_rep = rep(tool, length(common_names))
    tool_name_rep_total = append(tool_name_rep_total, tool_name_rep)  
}

# Matrix for celtype corr
celltype_proportions = data.frame(ref = vector1_common_all,
                                tool = vector2_common_all,
                                tool_names = tool_name_rep_total,
                                celltype_names = common_names_all,
                                tools_coor = tools_coor)

# Ordering 
celltype_proportions$tool_names = factor(celltype_proportions$tool_names , levels = list_of_tools)

celltype_proportions_vs = ggplot(data=celltype_proportions, aes(x=ref, y=tool, tool_name = tool_names, color=celltype_names))+ 
                    geom_point(size = 2) +
                    xlim(0, max_vector=max(c(celltype_proportions$ref, celltype_proportions$tool))) +
                    ylim(0, max_vector=max(c(celltype_proportions$ref, celltype_proportions$tool))) +
                    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
                    facet_wrap(~tool_names, ncol = length(list_of_tools)) +
                    labs(title = "Cell type proportion comparison", x = "Proportion in ref", y = "Proportion in tools") +
                    coord_fixed(ratio = 1) +
                    theme(plot.title = element_text(hjust = 0.5)) +
                    geom_text(data = celltype_proportions, aes(x = Inf, y = Inf, 
                                label = paste("Corr:", round(tools_coor, 2))), 
                                hjust = 1.1, vjust = 1.5, size = 2.5, color = "black", inherit.aes = FALSE) +
                    theme(legend.position = "right") +
                    guides(colour = guide_legend(ncol = 1, override.aes = list(size=4))) +
                    theme(
                        legend.key.size = unit(0.5, "lines"),
                        legend.text = element_text(size = 8)
                    )

pdf(paste0(output_path, "/images/", "Proportions_tool_celltype_vs_ref.pdf"), width = 10, height = 10)
print(celltype_proportions_vs)
dev.off()

######
# EXPRESSION PURITY
# 5th 
# Positive purity score for each tool compared to original data
positive_purity_tool_all = c()
positive_purity_ref_all = c()
celltype_names_all = c()
tool_name_rep_total = c()

for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")
    ref_obj_name = paste0(tool, "_ref_data")

    # Get positive purity values grouped by celltype
    positive_purity_score_tool = data.frame(colData(get(obj_name))) %>%
                                group_by(mean_celltype_correlation) %>%
                                summarize(average_value = mean(positive_F1, na.rm = TRUE)) 
    positive_purity_score_ref = data.frame(colData(spe)) %>%
                                group_by(mean_celltype_correlation) %>%
                                summarize(average_value = mean(positive_F1, na.rm = TRUE)) 

    # Common cells in both tool and ref
    common_names = intersect(positive_purity_score_tool$mean_celltype_correlation , positive_purity_score_ref$mean_celltype_correlation)
    
    # Subset base on common names 
    positive_purity_tool = positive_purity_score_tool[positive_purity_score_tool$mean_celltype_correlation %in% common_names,]$average_value
    positive_purity_ref = positive_purity_score_ref[positive_purity_score_ref$mean_celltype_correlation %in% common_names,]$average_value
    celltype_names = common_names

    # Agregate all tools values into one vector
    positive_purity_tool_all = append(positive_purity_tool_all, positive_purity_tool)
    positive_purity_ref_all = append(positive_purity_ref_all, positive_purity_ref)
    celltype_names_all = append(celltype_names_all, celltype_names)

    # Tool name rep 
    tool_name_rep = rep(tool, length(common_names))
    tool_name_rep_total = append(tool_name_rep_total, tool_name_rep)  
}

# Matrix for celtype corr
positive_purity_score_matrix = data.frame(ref = positive_purity_ref_all,
                                          tool = positive_purity_tool_all,
                                          celltype_names = celltype_names_all,
                                          tool_names = tool_name_rep_total)

# Ordering 
positive_purity_score_matrix$tool_names = factor(positive_purity_score_matrix$tool_names , levels = list_of_tools)

positive_purity_score_vs = ggplot(data=positive_purity_score_matrix, aes(x=ref, y=tool, tool_name = tool_names, color=celltype_names))+ 
                                    geom_point(size = 2) +
                                    xlim(0, max_vector=max(c(positive_purity_score_matrix$ref, positive_purity_score_matrix$tool))) +
                                    ylim(0, max_vector=max(c(positive_purity_score_matrix$ref, positive_purity_score_matrix$tool))) +
                                    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
                                    facet_wrap(~tool_names, ncol = length(list_of_tools)) +
                                    labs(title = "Positive purity score comparison", x = "Positive purity score (ref)", y = "Positive purity score tools") +
                                    coord_fixed(ratio = 1) +
                                    theme(plot.title = element_text(hjust = 0.5)) +
                                    theme(legend.position = "right") +
                                    guides(colour = guide_legend(ncol = 1, override.aes = list(size=4))) +
                                    theme(
                                        legend.key.size = unit(0.5, "lines"),
                                        legend.text = element_text(size = 8),
                                        axis.text.x = element_text(size = 7),
                                        axis.text.y = element_text(size = 7)   
                                    )

pdf(paste0(output_path, "/images/", "Positive_purity_tools_vs_nucleus_ref.pdf"), width = 10, height = 10)
print(positive_purity_score_vs)
dev.off()

######
# EXPRESSION PURITY
# 6th 
# Purity F1 score against average total transcripts per cell
list_of_tools_plus_ref = c("spe" , list_of_tools)
purity_score_all = c()
tool_total_transcripts_all = c()
tool_all = c()
for (tool in list_of_tools_plus_ref){
    print(tool)
    if (tool != "spe"){obj_name  = paste0(tool, "_spe")} else {obj_name = paste0(tool)}
    # Purity score computation for tool
    scaled_positive_F1_mean = mean(colData(get(obj_name))$positive_F1)
    scaled_negative_F1_mean = mean(colData(get(obj_name))$negative_F1)

    scaled_positive_F1_mean = (scaled_positive_F1_mean - min(colData(get(obj_name))$positive_F1)) / (max(colData(get(obj_name))$positive_F1) - min(colData(get(obj_name))$positive_F1))
    scaled_negative_F1_mean = (scaled_negative_F1_mean - min(colData(get(obj_name))$negative_F1)) / (max(colData(get(obj_name))$negative_F1) - min(colData(get(obj_name))$negative_F1))
    purity_score = 2 * (((1 - scaled_negative_F1_mean) * scaled_positive_F1_mean) / (1 - scaled_negative_F1_mean + scaled_positive_F1_mean))

    # Transcripts for tool
    tool_total_transcripts = mean(as.numeric(colData(get(obj_name))$total_transcripts))

    # Agregate all tools values into one vector
    purity_score_all = append(purity_score_all, purity_score)
    tool_total_transcripts_all = append(tool_total_transcripts_all, tool_total_transcripts)
    tool_all = append(tool_all, tool)
}

# Matrix for celtype corr
purity_score_matrix = data.frame(purity_score = purity_score_all,
                                total_transcripts = tool_total_transcripts_all,
                                tool_names = tool_all)

# Renaming and ordering 
purity_score_matrix$tool_names[purity_score_matrix$tool_names == "spe"] = "ref"
list_of_tools_plus_ref[list_of_tools_plus_ref == "spe"] = "ref"
purity_score_matrix$tool_names = factor(purity_score_matrix$tool_names , levels = list_of_tools_plus_ref)

purity_score_matrix = data.frame(purity_score_matrix)

purity_score_tools = ggplot(data=purity_score_matrix, aes(x=total_transcripts, y=purity_score, color=tool_names))+ 
                    geom_point(size = 8) +
                    labs(title = "", x = "Average total transcripts per cell", y = "Purity F1") +
                    theme(plot.title = element_text(hjust = 0.5)) +
                    geom_text(aes(label = tool_names), nudge_x = 8, nudge_y = max(purity_score_matrix$purity_score)/100) + 
                    theme(legend.position = "none")

pdf(paste0(output_path, "/images/", "Purity_score_tools_comparison.pdf"))
print(purity_score_tools)
dev.off()

######
# SPATIAL CHARACTERISTICS DIVERSITY
# 7th and 8th
# Enthropy celltype composition divertity
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")

    # 7th
    temp_tool_spatialCords = data.frame(get(obj_name)@metadata$CellSegOutput)

    # Compute centroids
    centroids = temp_tool_spatialCords %>%
    group_by(cell_id) %>%
    summarize(
        centroid_x = mean(coord_x),
        centroid_y = mean(coord_y)
    )
    temp_tool_spatialCords = cbind(centroids$centroid_x, centroids$centroid_y)
    colnames(temp_tool_spatialCords) = c("coord_x", "coord_y")

    # Add centroids to tools object
    assign_expr =  paste0(obj_name, " = get(obj_name)[,centroids$cell_id]")
    eval(parse(text = assign_expr))
    temp_env = get(obj_name)
    spatialCoords(temp_env) = temp_tool_spatialCords
    assign(obj_name, temp_env) 

    # Compute diversity over centroids
    assign(obj_name, calSpatialMetricsDiversity(get(obj_name), celltype = "mean_celltype_correlation", bin_width = 500)) 

    # Get into a matrix the diversity metrics
    diversity_df = get(obj_name)@metadata$CellSPA$spatialMetricsDiversity$results

    # Diversity matrix long version
    df_long = diversity_df %>%
    pivot_longer(cols = starts_with("cellTypeProp_"), 
                names_to = "cell_type", 
                values_to = "proportion")

    df_complete = df_long %>%
                    complete(x_bin = unique(df_long$x_bin), y_bin = unique(df_long$y_bin), fill = list(entropy = NA))

    # Generate the Spatial Regions plot with pie charts
    spatial_regions = ggplot(df_complete, aes(x = factor(1), y = proportion, fill = cell_type)) +
                            geom_bar(stat = "identity", width = 1) +
                            coord_polar(theta = "y") +
                            facet_wrap(~ interaction(x_bin, y_bin), nrow = length(unique(df_complete$y_bin)), ncol = length(unique(df_complete$x_bin))) +
                            theme_void() +
                            scale_fill_viridis_d() + # Change the scale to your preference
                            theme(aspect.ratio = 1, 
                                    panel.spacing = unit(-0.85, "lines"),  # Adjust panel spacing to control pie size
                                    strip.text = element_blank(),
                                    strip.background = element_blank(),
                                    legend.position = "none") + 
                            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))  # Remove extra margins

    pdf(paste0(output_path, "/images/", "Spatial_regions_", tool,".pdf"))
    print(spatial_regions)
    dev.off()

    # 8th
    # Reverse order to match both plots
    df_complete = df_complete %>%
    mutate(
        y_bin = factor(y_bin, levels = rev(levels(factor(y_bin))))
    )

    # Generate the Cell type composition diversity plot
    celltype_composition_diversity = ggplot(df_complete, aes(x = x_bin, y = y_bin, fill = entropy)) +
                                            geom_tile() +
                                            scale_fill_viridis_c() +
                                            theme(aspect.ratio = 1, 
                                                axis.text.x = element_blank(),
                                                axis.text.y = element_blank(),
                                                axis.ticks = element_blank())

    pdf(paste0(output_path, "/images/", "Composition_diversity_", tool, ".pdf"), width = 10, height = 10)
    print(celltype_composition_diversity)
    dev.off()
}

######
# SPATIAL CHARACTERISTICS DIVERSITY
# 9th
# Enthropy divertity
entropy_all = c()
cv_total_trnacripts_all = c()
tool_cor_all = c()
tool_name_rep_total = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")

    # 7th
    temp_tool_spatialCords = data.frame(get(obj_name)@metadata$CellSegOutput)

    # Compute centroids
    centroids = temp_tool_spatialCords %>%
    group_by(cell_id) %>%
    summarize(
        centroid_x = mean(coord_x),
        centroid_y = mean(coord_y)
    )
    temp_tool_spatialCords = cbind(centroids$centroid_x, centroids$centroid_y)
    colnames(temp_tool_spatialCords) = c("coord_x", "coord_y")

    # Subset spatialexperiment object for the current cells
    assign_expr = paste0(obj_name, "= get(obj_name)[, colnames(get(obj_name)) %in% centroids$cell_id]")
    eval(parse(text = assign_expr))

    # Add centroids to tools object
    temp_env = get(obj_name)
    spatialCoords(temp_env) = temp_tool_spatialCords
    assign(obj_name, temp_env) 

    # Compute diversity over centroids
    assign(obj_name, calSpatialMetricsDiversity(get(obj_name), celltype = "mean_celltype_correlation", bin_width = 500)) 

    # Get into a matrix the diversity metrics
    diversity_df = get(obj_name)@metadata$CellSPA$spatialMetricsDiversity$results
    
    # Extract values of interest
    entropy = diversity_df$entropy
    cv_total_trnacripts = diversity_df$cv_total_transcripts

    # Compute correlations to be added to the plots
    diversity_df$cv_total_transcripts[is.na(diversity_df$cv_total_transcripts)] = 0
    diversity_df$entropy[is.na(diversity_df$entropy)] = 0
    tool_cor = cor(diversity_df$entropy, diversity_df$cv_total_transcripts)
    tool_cor = rep(tool_cor, nrow(diversity_df))
    tool_cor_all = append(tool_cor_all, tool_cor)


    # Agregate all tools values into one vector
    entropy_all = append(entropy_all, entropy)
    cv_total_trnacripts_all = append(cv_total_trnacripts_all, cv_total_trnacripts)

    # Tool name rep 
    tool_name_rep = rep(tool, nrow(diversity_df))
    tool_name_rep_total = append(tool_name_rep_total, tool_name_rep)  
}

# Matrix for plotting
celltype_entropy= data.frame(entropy = entropy_all,
                                  cv_total_transcripts = cv_total_trnacripts_all,
                                  Tool_names = tool_name_rep_total,
                                  tool_cor_all = tool_cor_all)

# Ordering 
celltype_entropy$Tool_names = factor(celltype_entropy$Tool_names , levels = list_of_tools)


# Generate the plot
spatial_variation_2 = ggplot(celltype_entropy, aes(x = entropy, y = cv_total_transcripts, color = Tool_names)) +
                            geom_point() +
                            facet_wrap(~Tool_names, ncol = length(list_of_tools)) +
                            labs(title = "", x = "Entropy", y = "CV of total transcripts") +
                            geom_text(data = celltype_entropy, aes(x = Inf, y = Inf, 
                                label = paste("Corr:", round(tool_cor_all, 2))), 
                                hjust = 1.1, vjust = 1.5, size = 2.5, color = "black", inherit.aes = FALSE) +
                            theme(legend.position = "right") +
                            theme(aspect.ratio = 1)

pdf(paste0(output_path, "/images/", "Spatial_diversity_entropy.pdf"), width = 10, height = 10)
print(spatial_variation_2)
dev.off()

# 10th
# CV of elongation divertity (also add coords to spatialCoords slot of SE tools objects)
cv_of_elongation_all = c()
prop_celltype_all = c()
tool_cor_all2 = c()
tool_name_rep_total = c()
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")

    temp_tool_spatialCords = data.frame(get(obj_name)@metadata$CellSegOutput)

    # Compute centroids
    centroids = temp_tool_spatialCords %>%
    group_by(cell_id) %>%
    summarize(
        centroid_x = mean(coord_x),
        centroid_y = mean(coord_y)
    )
    temp_tool_spatialCords = cbind(centroids$centroid_x, centroids$centroid_y)
    colnames(temp_tool_spatialCords) = c("coord_x", "coord_y")

    # Add centroids to tools object
    temp_env = get(obj_name)
    spatialCoords(temp_env) = temp_tool_spatialCords
    assign(obj_name, temp_env) 

    # Compute diversity over centroids
    assign(obj_name, calSpatialMetricsDiversity(get(obj_name), celltype = "mean_celltype_correlation", bin_width = 500)) 

    # Get into a matrix the diversity metrics
    diversity_df = get(obj_name)@metadata$CellSPA$spatialMetricsDiversity$results
    
    # Extract values of interest
    cv_elongation = diversity_df$cv_elongation

    # Identify columns starting with "cellTypeProp_"
    cols = grep("^cellTypeProp_", colnames(diversity_df), value = TRUE)
    # Calculate the mean for each identified column
    means = sapply(diversity_df[, cols, drop = FALSE], mean, na.rm = TRUE)
    # Find the column with the highest mean
    highest_mean_col = cols[which.max(means)]

    prop_celltype = unlist(diversity_df[highest_mean_col] * 100)

    #Compute correlations to be added to the plots
    diversity_df[highest_mean_col][is.na(unlist(diversity_df[highest_mean_col]))] = 0
    diversity_df$cv_elongation[is.na(diversity_df$cv_elongation)] = 0

    tool_cor = cor(diversity_df$cv_elongation, diversity_df[highest_mean_col])
    tool_cor = rep(tool_cor, nrow(diversity_df))
    tool_cor_all2 = append(tool_cor_all2, tool_cor)

    # gregate all tools values into one vector
    cv_of_elongation_all = append(cv_of_elongation_all, cv_elongation)
    prop_celltype_all = append(prop_celltype_all, prop_celltype)

    # Tool name rep 
    tool_name_rep = rep(tool, nrow(diversity_df))
    tool_name_rep_total = append(tool_name_rep_total, tool_name_rep)  
}

# Matrix for plotting
celltype_elongation= data.frame(cv_of_elongation = cv_of_elongation_all,
                                  prop_celltype = prop_celltype_all,
                                  Tool_names = tool_name_rep_total,
                                  tool_cor_all2 = tool_cor_all2)

# Ordering 
celltype_elongation$Tool_names = factor(celltype_elongation$Tool_names , levels = list_of_tools)

# Generate the plot
spatial_variation_2 = ggplot(celltype_elongation, aes(x = prop_celltype, y = cv_of_elongation, color = Tool_names)) +
                            geom_point() +
                            facet_wrap(~Tool_names, ncol = length(list_of_tools)) +
                            labs(title = "", x = paste0(highest_mean_col, " %"), y = "CV of elongation") +
                            geom_text(data = celltype_elongation, aes(x = Inf, y = Inf, 
                                label = paste("Corr:", round(tool_cor_all2, 2))), 
                                hjust = 1.1, vjust = 1.5, size = 2.5, color = "black", inherit.aes = FALSE) +
                            theme(legend.position = "right") +
                            theme(aspect.ratio = 1)

pdf(paste0(output_path, "/images/", "Spatial_diversity_cv_elongation.pdf"), width = 10, height = 10)
print(spatial_variation_2)
dev.off()

######
# NEIGHBOURING CONTAMINATION
# 11th
# Neighbour purity

##################################################################################
# Add function code (list atributes changed for evade error from package fucntion)
calSpatialMetricsDiversity = function (spe, celltype, bin_width = 500) {
    if ("character" %in% class(celltype) & length(celltype) == 
        1) {
        celltype = colData(spe)[, celltype]
    }
    coord = SpatialExperiment::spatialCoords(spe)
    coord_x = coord[, 1]
    coord_y = coord[, 2]
    num_bin_x = ceiling(diff(range(coord_x))/bin_width)
    x_bin_groups = cut(coord_x, breaks = max(num_bin_x, 2))
    num_bin_y = ceiling(diff(range(coord_y))/bin_width)
    y_bin_groups = cut(coord_y, breaks = max(num_bin_y, 2))
    bin_groups = paste(x_bin_groups, y_bin_groups, sep = "|")
    df_bin = data.frame(x_bin_groups, y_bin_groups, bin_groups)
    colData(spe)[, colnames(df_bin)] = data.frame(x_bin_groups, 
        y_bin_groups, bin_groups)
    res = table(celltype, bin_groups)
    bin_entropy = apply(res, 2, function(x) entropy::entropy(x))
    df_entropy_bin = data.frame(bin = names(bin_entropy), entropy = bin_entropy)
    df_entropy_bin = cbind(df_entropy_bin, do.call(rbind, strsplit(as.character(df_entropy_bin$bin), 
        "\\|")))
    colnames(df_entropy_bin) = c("bin", "entropy", "x_bin", 
        "y_bin")
    df_entropy_bin$x_bin = factor(df_entropy_bin$x_bin, levels = levels(x_bin_groups))
    df_entropy_bin$y_bin = factor(df_entropy_bin$y_bin, levels = levels(y_bin_groups))
    celltype_prop = apply(res, 2, function(x) x/sum(x))
    celltype_prop = t(celltype_prop)
    colnames(celltype_prop) = paste("cellTypeProp", colnames(celltype_prop), 
        sep = "_")
    df_entropy_bin = cbind(df_entropy_bin, celltype_prop[as.character(df_entropy_bin$bin), 
        ])
    common_cell_metrics = c("total_transcripts", "total_genes", "cell_area") #This line changed from package
    for (i in common_cell_metrics) {
        if (i %in% c("total_transcripts", "total_genes", "cell_area")) {
            df_bin_stats = stats::aggregate(log10(colData(spe)[, 
                i]), list(bin_groups), function(x) sd(x)/mean(x))
            colnames(df_bin_stats) = c("bin", paste("cv", i, 
                sep = "_"))
            df_entropy_bin = merge(df_entropy_bin, df_bin_stats, 
                by = "bin")
            df_bin_stats = stats::aggregate((colData(spe)[, 
                i]), list(bin_groups), function(x) sd(x))
            colnames(df_bin_stats) = c("bin", paste("sd", i, 
                sep = "_"))
            df_entropy_bin = merge(df_entropy_bin, df_bin_stats, 
                by = "bin")
        }
        else {
            df_bin_stats = stats::aggregate((colData(spe)[, 
                i]), list(bin_groups), function(x) sd(x)/mean(x))
            colnames(df_bin_stats) = c("bin", paste("cv", i, 
                sep = "_"))
            df_entropy_bin = merge(df_entropy_bin, df_bin_stats, 
                by = "bin")
            df_bin_stats = stats::aggregate((colData(spe)[, 
                i]), list(bin_groups), function(x) sd(x))
            colnames(df_bin_stats) = c("bin", paste("sd", i, 
                sep = "_"))
            df_entropy_bin = merge(df_entropy_bin, df_bin_stats, 
                by = "bin")
        }
    }
    cor_cv = cor(df_entropy_bin[, "entropy"], df_entropy_bin[, 
        paste("cv", common_cell_metrics, sep = "_")], use = "complete.obs")[1, 
        ]
    cor_sd = cor(df_entropy_bin[, "entropy"], df_entropy_bin[, 
        paste("sd", common_cell_metrics, sep = "_")], use = "complete.obs")[1, 
        ]
    cor_res = c(cor_cv, cor_sd)
    spe@metadata$CellSPA$spatialMetricsDiversity = list(statistics = cor_res, 
        results = df_entropy_bin)
    spe = .add_metrics(spe, names(spe@metadata$CellSPA$spatialMetricsDiversity$statistics), 
        "bin_level")
    return(spe)
}

##################################################################################

list_of_tools_plus_ref = c("spe", list_of_tools)
celltype_1_distance_all = c()
celltype_2_distance_all = c()
tool_name_rep_total = c()
for (tool in list_of_tools_plus_ref){
    print(tool)
    if (tool != "spe"){obj_name  = paste0(tool, "_spe")} else {obj_name = paste0(tool)}

    # Add cellnames to rownames of spatial coords from tools objects
    temp_env = get(obj_name)
    rownames(spatialCoords(temp_env)) = colnames(get(obj_name))
    assign(obj_name, temp_env) 

    # Compute diversity over centroids
    assign(obj_name, calSpatialMetricsDiversity(get(obj_name), celltype = "mean_celltype_correlation", bin_width = 500)) 

    # Get into a matrix the diversity metrics
    diversity_df = get(obj_name)@metadata$CellSPA$spatialMetricsDiversity$results

    # Select cell types to study and get the negative markers
    one = "lung macrophage"
    second = "erythrocyte"
    nn_celltype_pair = c(one, second)

    # Select the genes that are known negative markers base on literature
    neg_markers = list(one = c("EPCAM"),
                        second = c("PTPRC"))
    neg_markers = setNames(neg_markers, nn_celltype_pair)

    # Compute distance from the nearest cell that contains selected gene markers, where the studied celltype selected are grouped by distance ranges             
    aaa = calNegMarkerVsDist(get(obj_name),
                            "mean_celltype_correlation",
                            nn_celltype_pair,
                            neg_markers,
                            distance_breaks = c(seq(0, 600, 50)))

    # Extract data
    neg_marker_data = aaa@metadata$CellSPA$negMarkerExprs_vs_dist

    # Convert to data frames
    celltype_1_distance_df = as.data.frame(neg_marker_data[1])
    celltype_2_distance_df = as.data.frame(neg_marker_data[2])

    # Select the range to visualize base on the data
    #distance_ranges = diversity_df$bin
    #upper_bounds = as.numeric(gsub("[\\(\\[]|,.*", "", distance_ranges))
    #max_value = max(upper_bounds)
    #range_selected = max_value/2

    # Add distance column
    celltype_1_distance_df$distance = rownames(celltype_1_distance_df)
    celltype_2_distance_df$distance = rownames(celltype_2_distance_df)

    # Reshape the data for ggplot2
    celltype_1_distance = tidyr::gather(celltype_1_distance_df, key = "gene", value = "expression", -distance)
    celltype_2_distance = tidyr::gather(celltype_2_distance_df, key = "gene", value = "expression", -distance)

    celltype_1_distance$distance = factor(celltype_1_distance$distance, levels = celltype_1_distance$distance)   
    celltype_2_distance$distance = factor(celltype_2_distance$distance, levels = celltype_2_distance$distance)  

    # Select first 10 distances
    celltype_1_distance = celltype_1_distance[1:5,]
    celltype_2_distance = celltype_2_distance[1:5,]

    # Add tool name
    tool_name_rep = rep(tool, nrow(celltype_1_distance))
    tool_name_rep_total = append(tool_name_rep_total, tool_name_rep)  

    celltype_1_distance_all = rbind(celltype_1_distance_all, celltype_1_distance)     
    celltype_2_distance_all = rbind(celltype_2_distance_all, celltype_2_distance)   
}

# Renaming 
tool_name_rep_total[tool_name_rep_total == "spe"] = "ref"
list_of_tools_plus_ref[list_of_tools_plus_ref == "spe"] = "ref"

# Add tools name vector to tools distance matrix
celltype_1_distance_all$tool_names = tool_name_rep_total
celltype_2_distance_all$tool_names = tool_name_rep_total

# Ordering 
celltype_1_distance_all$tool_names = factor(celltype_1_distance_all$tool_names , levels = list_of_tools_plus_ref)
celltype_2_distance_all$tool_names = factor(celltype_2_distance_all$tool_names , levels = list_of_tools_plus_ref)

neighbour_purity_1= ggplot(celltype_1_distance_all, aes(x = distance, y = expression, color = tool_names, group = interaction(gene, tool_names))) +
  geom_line() +
  geom_point() +
  labs(title = paste0(unique(celltype_1_distance_all$gene) ," in ", nn_celltype_pair[1]), x = "Euclidean distance", y = "% cells expressing negative marker") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

pdf(paste0(output_path, "/images/", "Neighbour_purity_",one,".pdf"), width = 10, height = 10)
print(neighbour_purity_1)
dev.off()

neighbour_purity_2 = ggplot(celltype_2_distance_all, aes(x = distance, y = expression, color = tool_names, group = interaction(gene, tool_names))) +
  geom_line() +
  geom_point() +
  labs(title = paste0(unique(celltype_2_distance_all$gene) , " in ", nn_celltype_pair[2]), x = "Euclidean distance", y = "% cells expressing negative marker") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

pdf(paste0(output_path, "/images/", "Neighbour_purity_",second,".pdf"), width = 10, height = 10)
print(neighbour_purity_2)
dev.off()



######
# UMAP VISUALIZATION
# 12th
# UMAPS of each tool base on spatial rna information
list_of_tools_plus_ref = c("spe", list_of_tools)
umap_coords_all = c()
celltype_all = c()
tool_name_rep_total = c()
for (tool in list_of_tools_plus_ref){
    print(tool)
    if (tool != "spe"){obj_name  = paste0(tool, "_spe")} else {obj_name = paste0(tool)}

    # Extract umap dim for each tool
    umap_coords = reducedDim(get(obj_name), "UMAP")
    umap_coords_all = rbind(umap_coords_all, umap_coords)

    # Add celltype for each tool
    celltype = get(obj_name)$mean_celltype_correlation
    celltype_all = append(celltype_all, celltype)

    # Add tool name
    tool_name_rep = rep(tool, nrow(umap_coords))
    tool_name_rep_total = append(tool_name_rep_total, tool_name_rep) 
}

# Renaming 
tool_name_rep_total[tool_name_rep_total == "spe"] = "ref"
list_of_tools_plus_ref[list_of_tools_plus_ref == "spe"] = "ref"

# Matrix for plotting
umap_dim_matrix = data.frame(umap_dims_1 = umap_coords_all[,1],
                             umap_dims_2 = umap_coords_all[,2],
                             tool_names = tool_name_rep_total,
                             celltypes = celltype_all)

# Ordering 
umap_dim_matrix$tool_names = factor(umap_dim_matrix$tool_names , levels = list_of_tools_plus_ref)
umap_dim_matrix$tool_names = factor(umap_dim_matrix$tool_names , levels = list_of_tools_plus_ref)

umap_plot = ggplot(umap_dim_matrix, aes(x = umap_dims_1, y = umap_dims_2)) +
    geom_point(aes(color = celltypes), size = 1) +
    facet_wrap(~tool_names, ncol = 2, nrow = 4) +
    labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
    theme(aspect.ratio = 1, legend.position = "right")  +
    guides(color = guide_legend(ncol = 1))   


# Display the plot
pdf(paste0(output_path, "/images/", "UMAP_tools.pdf"), width = 10, height = 15)
print(umap_plot)
dev.off()

######
# PATCH VISUALIZATION WITH SEGMENTATIONS
#1 3th
# Plot images with segmentation masks
plot_image_and_mask = function(image_path, random_patch, x_relocation, y_relocation) {
    # Read the image
    img = image_read(image_path)

    # Read the mask data
    mask_data = random_patch
    colnames(mask_data) = c("x","y","cell_id")

    # Convert the image to a ggplot-friendly format
    img_raster = as.raster(img)
    img_data = data.frame(
        x = rep(1:ncol(img_raster), each = nrow(img_raster)),
        y = rep(1:nrow(img_raster), times = ncol(img_raster)),
        fill = as.vector(img_raster)
    )

    # To print the images in the correct image position
    img_data$x = max(img_data$x) - img_data$x
    mask_data$x = mask_data$x - x_relocation
    mask_data$y = mask_data$y - y_relocation
    mask_data$y = 3600 - mask_data$y

    # Aggregate points into polygons
    mask_data_grouped = mask_data %>%
        group_by(cell_id) %>%
        summarize(x = list(x), y = list(y)) %>%
        ungroup()
    
    # Create convex hulls around each group
    convex_hulls = mask_data_grouped %>%
        rowwise() %>%
        do({
            df <- data.frame(x = unlist(.$x), y = unlist(.$y))
            hull <- chull(df$x, df$y)
            df <- df[hull, ]
            df$cell_id <- .$cell_id
            df
        }) %>%
        ungroup()

    # # Plot the image
    # ggplot() +
    #     geom_raster(data = img_data, aes(x = y, y = x, fill = fill)) +
    #     scale_fill_identity() +
    #     theme_void() +
    #     coord_fixed() +
    #     geom_point(data = mask_data, aes(x = x, y = y, color = factor(cell_id)),alpha = 0.05) +
    #     scale_color_viridis_d() + # Customize colors as needed
    #     theme(legend.position = "none")

    # Plot the image
    ggplot() +
        geom_raster(data = img_data, aes(x = y, y = x, fill = fill), show.legend = FALSE) + # Remove legend for image
        scale_fill_identity() +
        theme_void() +
        coord_fixed() +
        geom_polygon(data = convex_hulls, aes(x = x, y = y, color = factor(cell_id)), alpha = 0.5) +
        scale_color_viridis_d() + 
        theme(legend.position = "none")
}

count = 0
for (tool in list_of_tools){
    print(tool)
    obj_name  = paste0(tool, "_spe")

    # Get all patches for each tool
    output_paths = process_paths(tool_used = tool)

    # Extract random transcriptomic patch information
    extracted_numbers_1 = gsub(".*\\/([0-9:]+[0-9:]+)\\.txt$", "\\1", output_paths)
    extracted_numbers_2 = gsub(".*\\/([0-9]+:[0-9]+_[0-9]+:[0-9]+)\\.txt$", "\\1", output_paths)
    extracted_numbers_3 = gsub(".*spot2nucl_([0-9:]+)\\.txt$", "\\1", output_paths)
    extracted_numbers_4 = gsub(".*spot2cell_([0-9:]+)\\.txt$", "\\1", output_paths)

    if (tool == "bidcell"){
        # Iterate through each file in the initial directory
        extracted_numbers_1 = c()
        for (file in output_paths) {
            # List files in another subdirectory
            another_subdirectory_files <- list.files(file, full.names = TRUE)[8]
            #
            # List specific files with a given pattern or index
            output_paths_n = list.files(another_subdirectory_files, full.names = TRUE)
            output_paths_n = list.files(output_paths_n, full.names = TRUE)[3]
            output_paths_n_tiff = list.files(output_paths_n, full.names = TRUE)[3]
            #
            # Add these files to the list of all file paths
            extracted_numbers_1 <- append(extracted_numbers_1, output_paths_n_tiff)
        }
    }

    # Find the minimum length and return the corresponding string
    find_shortest <- function(...) {
        inputs <- list(...)
        lengths <- sapply(inputs, nchar)
        min_index <- which.min(lengths)
        return(inputs[[min_index]])
    }
    if (tool == "bidcell"){
        extracted_numbers <- mapply(function(...) {
            do.call(find_shortest, list(...))
        }, extracted_numbers_1)    
    } else {
        extracted_numbers <- mapply(function(...) {
            do.call(find_shortest, list(...))
        }, extracted_numbers_1, extracted_numbers_2, extracted_numbers_3,extracted_numbers_4)       
    }
    extracted_numbers <- unique(extracted_numbers)

    # Randomize only once (use the same one for all tools)
    count = count + 1
    if (count == 1){
        random_patch_selection = sample(extracted_numbers, 1)
        random_patch_selection = basename(random_patch_selection)
        random_patch_selection_split = as.numeric(unlist(strsplit(gsub("[^0-9]+", " ", random_patch_selection), "\\s+")[[1]]))
    }
    
    random_patch_transcripts = get(obj_name)@metadata$CellSegOutput[(get(obj_name)@metadata$CellSegOutput['coord_x'] >= random_patch_selection_split[2]+1) & 
                                (get(obj_name)@metadata$CellSegOutput['coord_x'] <= random_patch_selection_split[2]+random_patch_selection_split[3]-1) & 
                                (get(obj_name)@metadata$CellSegOutput['coord_y'] >= random_patch_selection_split[1]+1) & 
                                (get(obj_name)@metadata$CellSegOutput['coord_y'] <= random_patch_selection_split[1]+random_patch_selection_split[4]-1),]
    
    # Extract the corresponding image path image for random transcriptome selected
    random_patch_image_path = "/ibex/scratch/projects/c2169/Xabier/spatial_transcriptomics/datasets/breast_data/crop_images"
    png_files <- list.files(path = random_patch_image_path, pattern = "\\.png$", full.names = TRUE)
    pattern <- paste0(random_patch_selection_split[1], "[:_]", random_patch_selection_split[2], "[:_]", random_patch_selection_split[3], "[:_]", random_patch_selection_split[4])
    matching_png_files = png_files[grep(pattern, png_files)][1]     

    pdf(paste0(output_path, "/images/", "Segmentation_visualization_",tool,".pdf"), width = 10, height = 10)
    print(plot_image_and_mask(image_path = matching_png_files, random_patch = random_patch_transcripts, x_relocation = random_patch_selection_split[2], y_relocation = random_patch_selection_split[1]))
    dev.off()
}


######
# OVERALL RANKING FOR SPECIFIC DATASET
# Metrics 1st panel

# Function to scale scores to the range [0, 1]
scale_scores <- function(scores) {
  # Find the minimum and maximum scores
  min_score <- min(unlist(scores))
  max_score <- max(unlist(scores))
  
  # Scale the scores
  scaled_scores <- lapply(scores, function(score) {
    (score - min_score) / (max_score - min_score)
  })
  scaled_scores = unlist(scaled_scores)
  return(scaled_scores)
}

# Function to compute silhouette score
compute_silhouette_scores <- function(umap_data, celltypes_col, toolname_col, dimensions = c("umap_dims_1", "umap_dims_2")) {
  # List to store silhouette scores for each tool
  silhouette_scores <- list()
  
  # Get unique tools
  tools <- unique(umap_data[[toolname_col]])
  
  # Loop over each tool
  for (tool in tools) {
    # Filter data for the current tool
    tool_data <- umap_data %>% filter(!!sym(toolname_col) == tool)
    
    # Extract dimensions and cell types
    umap_dims <- tool_data[, dimensions]
    celltypes <- tool_data[[celltypes_col]]
    
    # Convert cell types to numeric labels for clustering
    cluster_labels <- as.numeric(factor(celltypes))
    
    # Compute the distance matrix
    dist_matrix <- dist(umap_dims)
    
    # Compute silhouette widths
    silhouette_scores_tool <- silhouette(cluster_labels, dist_matrix)
    
    # Average silhouette width
    avg_silhouette_score <- mean(silhouette_scores_tool[, "sil_width"])
    
    # Store the score in the list
    silhouette_scores[[tool]] <- avg_silhouette_score
  }
  
  return(silhouette_scores)
}



#############################
#############################
# Create empty matrix for ranking metrics
#############################
#############################

# Changed required to fit size
purity_score_matrix_filtered = purity_score_matrix[purity_score_matrix["tool_names"] != "ref",]
purity_score_matrix_filtered$tool_names <- droplevels(purity_score_matrix_filtered$tool_names)
#
umap_dim_matrixfiltered = umap_dim_matrix[umap_dim_matrix["tool_names"] != "ref",]
umap_dim_matrixfiltered$tool_names <- droplevels(umap_dim_matrixfiltered$tool_names)
# Subsampling is needed as the data could be too large to compute silhouette and rise C erro
umap_dim_matrixfiltered <- data.frame(umap_dim_matrixfiltered %>%
                                group_by(tool_names, celltypes) %>%
                                sample_frac(0.4))

######
# BUILD RANKING MATRIX
ranking_metrics <- data.frame(
    # 1st panel
    # Overall characteristics
    hq_segmented_cells = scale_scores(ncell_comp$cells),
    pct_transcripts_assigned = scale_scores(transcripts_pct_comp$pct),
    # Cell-lvl QC metrics
    total_transcripts_per_cell = scale_scores((transcripts_total_comp %>% group_by(tool) %>% summarize(total = mean(total)))$total),
    total_genes_per_cell = scale_scores((genes_total_comp %>% group_by(tool) %>% summarize(total = mean(total)))$total),
    # Gene-lvl QC metrics
    pct_assigned_genes = scale_scores((pct_assigned_genes %>% group_by(tool) %>% summarize(mean_ratio = mean((pct_assig_genes / ref_pct_assig_genes)[is.finite(pct_assig_genes / ref_pct_assig_genes)], na.rm = TRUE)))$mean_ratio),
    # Cell morphology
    cor_transcripts_vs_ref = c(NA,scale_scores(transcipts_total_comp$tool_correlations)),
    cor_ct_elongation_vs_transcripts = c(NA, scale_scores(transcipts_total_comp$tool_correlations * transcipts_total_comp$tool_average_transcripts)),
    # 2nd panel
    # Expression and cell composition similarity with reference
    avg_diag_exp_vs_ref = scale_scores(avg_diag_value_all),
    ct_correlations = scale_scores(celltype_correlations$sum_value * celltype_correlations$tools_average_total_transcripts),
    expression_cor_vs_total_transcripts = scale_scores(celltype_correlations2$sum_value * celltype_correlations2$tools_average_total_transcripts),
    cor_ct_proportion = scale_scores((celltype_proportions %>% group_by(tool_names) %>% summarize(total = mean(tools_coor)))$total),
    # Expression purity
    positive_purity_score = scale_scores((positive_purity_score_matrix %>% group_by(tool_names) %>% summarize(total = mean(tool / ref)))$total),
    purity_score_F1 = scale_scores((purity_score_matrix_filtered %>% group_by(tool_names) %>% summarize(total = mean(purity_score/total_transcripts)))$total),
    # Spatial characteristics diversity
    cor_cv_trans_vs_entropy = scale_scores((celltype_entropy %>% group_by(Tool_names) %>% summarize(total = mean(tool_cor_all)))$total),
    #cor_cv_elongation_vs_celltype = scale_scores((celltype_elongation %>% group_by(Tool_names) %>% summarize(total = mean(tool_cor_all)))$total),
    # UMAPS representation consistency with celltype 
    umap_silhouette_scores = scale_scores(compute_silhouette_scores(umap_dim_matrixfiltered, "celltypes", "tool_names"))
)

# Plot dataset overall ranking
dataset_overall_scores = rowMeans(ranking_metrics,na.rm = TRUE)
dataset_overall_scores_df <- data.frame(
    tool_names = names(dataset_overall_scores),
    score = as.numeric(dataset_overall_scores)
)

# Create a ggplot color palette to match facet_wrap used in other plots
default_colors = c(
  "#F8766D", # Red
  "#B79F00", # Gold
  "#00BA38", # Green
  "#00BFC4", # Teal
  "#619CFF", # Blue
  "#F564E3", # Purple
  "#FF61CC" # Pink
)

# Ordering 
dataset_overall_scores_df$tool_names = factor(dataset_overall_scores_df$tool_names , levels = list_of_tools)

# Barplot with tools overall ranking in specific dataset
barplot <- ggplot(dataset_overall_scores_df, aes(x = tool_names, y = score, fill = tool_names)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = default_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Score", title = NULL)

pdf(paste0(output_path, "/images/", "Overall_ranking_barplot.pdf"), width = 10, height = 10)
print(barplot)
dev.off()

# Compute correlation matrix
corr_ranking_metrics = cor(ranking_metrics, use = "pairwise.complete.obs")

# Reorder base on hierarchical clustering
hc = hclust(as.dist(1 - abs(corr_ranking_metrics)))
corr_ranking_metrics = corr_ranking_metrics[hc$order, hc$order]

# Transform the correlation matrix to long format
melted_corr_ranking_metrics = melt(corr_ranking_metrics)

# Create the heatmap using ggplot2
cor_ranking =  ggplot(melted_corr_ranking_metrics, aes(x = Var1, y = Var2, fill = value)) +
                    geom_tile() +
                    scale_fill_viridis(option = "M", name = "Correlation") + 
                    theme_minimal() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    labs(x = NULL, y = NULL, title = "Correlation Heatmap")

pdf(paste0(output_path, "/images/", "Overall_ranking_cor_heatmap.pdf"), width = 10, height = 10)
print(cor_ranking)
dev.off()

######
# BUILD TABLE OF METRIX PLOT

######
# BUILD RANKING MATRIX BY CATEGORIES
overall_characteristics = rowMeans(ranking_metrics[,c("hq_segmented_cells","pct_transcripts_assigned")])
cell_level_qc_metrics = rowMeans(ranking_metrics[,c("total_transcripts_per_cell","total_genes_per_cell")])
gene_level_pct_assigned_genes = ranking_metrics["pct_assigned_genes"]
cell_morphology_metrics = rowMeans(ranking_metrics[,c("cor_transcripts_vs_ref","cor_ct_elongation_vs_transcripts")])
exp_cell_composition_similarity = rowMeans(ranking_metrics[,c("avg_diag_exp_vs_ref","ct_correlations","expression_cor_vs_total_transcripts","cor_ct_proportion")])
exp_purity = rowMeans(ranking_metrics[,c("positive_purity_score","purity_score_F1")])
spatial_characteristics = ranking_metrics["cor_cv_trans_vs_entropy"]
umaps_silhouette = ranking_metrics["umap_silhouette_scores"]

ranking_metrics = data.frame(overall_characteristics,cell_level_qc_metrics,gene_level_pct_assigned_genes,cell_morphology_metrics,exp_cell_composition_similarity,exp_purity,spatial_characteristics,umaps_silhouette)

# Plot dataset overall ranking
dataset_overall_scores = rowMeans(ranking_metrics)
dataset_overall_scores_df <- data.frame(
    tool_names = names(dataset_overall_scores),
    score = as.numeric(dataset_overall_scores)
)

# Create a ggplot color palette to match facet_wrap used in other plots
default_colors = c(
  "#F8766D", # Red
  "#B79F00", # Gold
  "#00BA38", # Green
  "#00BFC4", # Teal
  "#619CFF", # Blue
  "#F564E3", # Purple
  "#FF61CC" # Pink
)

# Ordering 
dataset_overall_scores_df$tool_names = factor(dataset_overall_scores_df$tool_names , levels = list_of_tools)

# Barplot with tools overall ranking in specific dataset
barplot <- ggplot(dataset_overall_scores_df, aes(x = tool_names, y = score, fill = tool_names)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = default_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Score", title = NULL)

pdf(paste0(output_path, "/images/", "Overall_ranking_barplot_by_categories.pdf"), width = 10, height = 10)
print(barplot)
dev.off()

# Compute correlation matrix
corr_ranking_metrics = cor(ranking_metrics, use = "pairwise.complete.obs")

# Reorder base on hierarchical clustering
hc = hclust(as.dist(1 - abs(corr_ranking_metrics)))
corr_ranking_metrics = corr_ranking_metrics[hc$order, hc$order]

# Transform the correlation matrix to long format
melted_corr_ranking_metrics = melt(corr_ranking_metrics)

# Create the heatmap using ggplot2
cor_ranking =  ggplot(melted_corr_ranking_metrics, aes(x = Var1, y = Var2, fill = value)) +
                    geom_tile() +
                    scale_fill_viridis(option = "M", name = "Correlation") +  # Apply viridis color palette
                    theme_minimal() +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    labs(x = NULL, y = NULL, title = "Correlation Heatmap")

pdf(paste0(output_path, "/images/", "Overall_ranking_cor_heatmap_by_categories.pdf"), width = 10, height = 10)
print(cor_ranking)
dev.off()

# Radar plot with categories (we can see the wheight of each category)
# Add maximum and minimum rows (for radar chart scaling)
max_values = rep(1.0, ncol(ranking_metrics))
min_values = rep(0.0, ncol(ranking_metrics))
ranking_metrics = rbind(max_values, min_values, ranking_metrics)

# Create a ggplot color palette to match facet_wrap used in other plots
default_colors = c(
  "#F8766D", # Red
  "#B79F00", # Gold
  "#00BA38", # Green
  "#00BFC4", # Teal
  "#619CFF", # Blue
  "#F564E3", # Purple
  "#FF61CC" # Pink
)

# Plot radar
pdf(paste0(output_path, "/images/", "Overall_ranking_radar_by_categories.pdf"), width = 10, height = 10)
# Adjust margins
par(mar = c(5, 5, 5, 5), oma = c(0, 0, 0, 0), xpd = TRUE)
radarchart(ranking_metrics, axistype = 1,
           # Customize the appearance
           plty = 1, # Line type
           # Define axis labels, scaling
           cglcol = "grey", # Grid color
           cglty = 1, # Grid line type
           axislabcol = "black", # Axis label color
           caxislabels = seq(0, 1, 0.25), # Axis labels
           cglwd = 1, # Grid line width
           vlcex = 0.60, # Increase this value to move labels further outside
           pcol = default_colors[1:8], # Line colors for each entity
           pfcol = scales::alpha(default_colors[1:8], 0.1))
# Add a legend
legend(x = 1.1, y = 1, 
       legend = row.names(ranking_metrics)[-c(1, 2)], # Exclude "Max" and "Min"
       col = default_colors,
       lty = 1, # Line type
       lwd = 2, # Line width
       bty = "n", # No box around the legend
       cex = 1.0) # Text size
dev.off()

# Save results with image and rank table
save.image(paste0(output_path, "/last_plot.RData"))
