Nextflow pipeline for segmentation methods
===========

This folder contains the Nextflow baseline framework for running segmentation methods in a unified way.

User steps:
- **Conda environments**: Ensure you have required conda environments before running.
- **Nextflow directory**: Clone the Nextflow directory provided in this repository.
- **Change output directory**: Change the output directory where outcomes should be stored in main.nf file.
- **Add input data**: Add to the input_data folder the two required files, the trsanscripts and the image.
- **Run the Nextflow**: Run the command 
      ```
      $ nextflow run .../nextflow/main.nf
      ```

The pipeline will preprocess, generate the patches and run the added segmentation tools.
