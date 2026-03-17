# Spatial Transcriptomics 
## Pipeline
A Snakemake-based spatial transcriptomics workflow for processing and analyzing 10x Genomics **Visium** and **Xenium** data. The repository is organized to support preprocessing outputs from Space Ranger / Xenium-style inputs and downstream analysis in **Seurat**. At present, the active workflow centers on per-sample Seurat object generation and QC/visualization outputs

### Running the Pipeline
1. Clone this repository and move into the directory: 
```
git clone https://github.com/xtinacomo/card_spatialtranscriptomics.git
cd card_spatialtranscriptomics
```

### Inputs:
2. Edit inputs/samples.csv file to include sample information

![Screenshot 2025-04-15 at 1 25 26 PM](https://github.com/user-attachments/assets/8c57334e-b357-41e1-b73a-cab2d290b6f8)

3. Add your data to the data folder and each sample has its' own folder:

![Screenshot 2025-04-15 at 1 32 59 PM](https://github.com/user-attachments/assets/7166696d-0bc7-4284-8196-60eb9fe86827)

5. For Visium samples, customize the snakefile to include the path to your transcriptome file 
6. Run snakemake.sh file

```
bash snakemake.sh
```

### Outputs:

![Screenshot 2025-04-15 at 1 34 58 PM](https://github.com/user-attachments/assets/6afdeb81-9fd4-44a5-a69a-eaddd795effd)

## Repository structure

```text
card_spatialtranscriptomics/
├── Snakefile
├── config.yaml
├── snakemake.sh
├── README.md
├── inputs/
│   └── samples.csv
├── logs/
├── probesets/
└── scripts/
    ├── seurat_process.R
    └── merge_seurat.R

