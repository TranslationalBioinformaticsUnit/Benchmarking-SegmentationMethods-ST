#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Input and output parameters
params.output1 = "/ibex/scratch/projects/c2169/Xabier/nextflow/1_processed_files"
params.output2 = "/ibex/scratch/projects/c2169/Xabier/nextflow/2_image_patches"
params.output3 = "/ibex/scratch/projects/c2169/Xabier/nextflow/3_transcript_patches"
params.output4 = "/ibex/scratch/projects/c2169/Xabier/nextflow/4__SCS_segmented_cells"
params.output5 = "/ibex/scratch/projects/c2169/Xabier/nextflow/5_Cellpose_segmented_cells"

// Process parameters for script1
params.transcripts_input = "/ibex/scratch/projects/c2169/Xabier/nextflow/input_data/slice1_points.csv"
params.image_input = "/ibex/scratch/projects/c2169/Xabier/nextflow/input_data/1_image.tif"
params.output_path = "Lung"
params.color_code = "rgb"

// Process parameters for script2
params.patch_size = 3600
params.subset = "False"

// Process parameters for script5
params.cellpose_model = "cyto3"
params.cellpose_restore_type = "denoise_cyto3"

/////////////////////////////////////////////////////////////////////////////
// Process 1
process run_script1 {
    conda '/home/iribarxm/miniconda3/envs/scs'

    output:
    path "${params.output_path}/transcripts_formated_full.csv"
    path "${params.output_path}/transcripts_formated.csv"
    path "${params.output_path}/Image_formated.tiff"

    publishDir params.output1, mode: 'copy'  // Copy outputs to the output1 directory

    script:
    """
    python /ibex/scratch/projects/c2169/Xabier/nextflow/bin/script1.py \
        ${params.transcripts_input} \
        ${params.image_input} \
        ${params.output_path} \
        ${params.color_code}
    """
}

// Process 2
process run_script2 {
    conda '/home/iribarxm/miniconda3/envs/scs'

    input:
    path bin_file
    path tiff_image

    output:
    path "${params.output_path}/*.png"
    path "${params.output_path}/*.tif"
    path "${params.output_path}/full_image/*"

    publishDir params.output2, mode: 'copy'

    script:
    """
    python /ibex/scratch/projects/c2169/Xabier/nextflow/bin/script2.py \
        ${bin_file} \
        ${tiff_image} \
        ${params.patch_size} \
        ${params.output_path} \
        ${params.subset} \
        ${params.color_code}
    """
}

// Process 3
process run_script3 {
    conda '/home/iribarxm/miniconda3/envs/scs'

    input:
    path bin_file

    output:
    path "${params.output_path}/*.csv"

    publishDir params.output3, mode: 'copy'

    script:
    """
    python /ibex/scratch/projects/c2169/Xabier/nextflow/bin/script3.py \
        ${bin_file} \
        ${params.output_path} \
        ${params.patch_size}
    """
}

// Process 4
process run_script4 {
    conda '/home/iribarxm/miniconda3/envs/scs'

    input:
    path bin_file
    path tiff_image

    output:
    path "${params.output_path}/results/*"

    publishDir params.output4, mode: 'copy'

    script:
    """
    python /ibex/scratch/projects/c2169/Xabier/nextflow/bin/script4.py \
        ${bin_file} \
        ${tiff_image} \
        ${params.output_path}
    """
}

// Process 5
process run_script5 {
    conda '/home/iribarxm/miniconda3/envs/cellpose'

    input:
    path image_patches

    output:
    path "${params.output_path}/*.txt"
    path "${params.output_path}/*.pdf"

    publishDir params.output5, mode: 'copy'

    script:
    """
    python /ibex/scratch/projects/c2169/Xabier/nextflow/bin/script5.py \
        ${params.output2} \
        ${params.cellpose_model} \
        ${params.cellpose_restore_type} \
        ${params.output_path} 
    """
}

/////////////////////////////////////////////////////////////////////////////
// Workflow definition
workflow {
    // Process files
    output1_files = run_script1()

    // Generate image patches
    output2_files = run_script2(bin_file = output1_files[1], tiff_image = output1_files[2])

    // Generate transcript patches
    // output3_files = run_script3(bin_file = output1_files[0])

    // Run SCS segmentation
    // output4_files = run_script4(bin_file = output1_files[1], tiff_image = output1_files[2])

    // Run Cellpose segmentation
    output5_files = run_script5(output2_files[0])
}
