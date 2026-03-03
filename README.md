# Spatial Transcriptomics 
## Pipeline
This pipeline processes 10x Genomics Visium and Xenium spatial tranascriptomics data from fastq files using [Snakemake](https://snakemake.readthedocs.io), SpaceRanger and [Seurat](https://satijalab.org/seurat/).

### Requirements
- Snakemake ≥ v7
- R ≥ v4.2
- R packages: Seurat, ggplot2, dplyr, patchwork
- Spaceranger (for Visium samples)
- Linux environment with SLURM (Biowulf-friendly)

### Running the Pipeline
1. Clone this repository and move into the directory: 
```
git clone https://github.com/NIH-CARD/card-unified-workflow/SpatialTranscriptomics.git
cd SpatialTranscriptomics
```

### Inputs:
2. Edit inputs/samples.csv file to include sample information

![Screenshot 2025-04-15 at 1 25 26 PM](https://github.com/user-attachments/assets/8c57334e-b357-41e1-b73a-cab2d290b6f8)

3. Add your data to the data folder and each sample has its' own folder:

![Screenshot 2025-04-15 at 1 32 59 PM](https://github.com/user-attachments/assets/7166696d-0bc7-4284-8196-60eb9fe86827)

5. For Visium samples, customize the snakefile to include the path to your transcriptome file 
6. Run snakemake.sh file

### Outputs:

![Screenshot 2025-04-15 at 1 34 58 PM](https://github.com/user-attachments/assets/6afdeb81-9fd4-44a5-a69a-eaddd795effd)

