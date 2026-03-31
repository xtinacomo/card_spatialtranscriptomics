import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

# File locations
data_dir = 'data/'
#results_dir = 'results/'

# Load sample metadata
samples_df = pd.read_csv('inputs/samples.csv')

#Grab samples,slid_id, and capture_area from samples.csv
SAMPLES = samples_df["sample"].astype(str).tolist()
SLIDE_DICT = samples_df.set_index("sample")["slide_id"].to_dict()
AREA_DICT = samples_df.set_index("sample")["area"].to_dict()
PLATFORM_DICT = samples_df.set_index("sample")["platform"].to_dict()

print(f"SAMPLES: {SAMPLES}")
print(f"SLIDE_DICT: {SLIDE_DICT}")
print(f"AREA_DICT: {AREA_DICT}")
print(f"PLATFORM_DICT: {PLATFORM_DICT}")

"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Final targets for both platforms
rule all:
    input:
        expand("results/seurat/{sample}/seurat.rds", sample=SAMPLES),
        expand("results/seurat/{sample}/spatial_plot.pdf", sample=SAMPLES)

def get_seurat_input(wildcards):
    platform = PLATFORM_DICT[wildcards.sample].lower()
    if platform == "visium":
        return f"results/{wildcards.sample}/outs/filtered_feature_bc_matrix.h5"
    elif platform == "xenium":
        return f"data/{wildcards.sample}/cell_feature_matrix.h5"
    else:
        raise ValueError(f"Unknown platform '{platform}' for sample {wildcards.sample}")

# SpaceRanger count for Visium
rule spaceranger_count:
    input:
        fastqs = lambda wc: os.path.join(data_dir, f"{wc.sample}_fastqs"),
        image  = lambda wc: os.path.join(data_dir, f"{wc.sample}_fastqs", f"{wc.sample}_image.tif")
    output:
        "results/{sample}/outs/filtered_feature_bc_matrix.h5"
    params:
        slide = lambda wc: SLIDE_DICT[wc.sample],
        area  = lambda wc: AREA_DICT[wc.sample],
        transcriptome= config["transcriptome"],
        create_bam   = config["create_bam"],
        nucleus_segmentation = lambda wc: f"--nucleus-segmentation={config['nucleus_segmentation']}" if config.get("nucleus_segmentation") not in [False, "false", "False", None] else "",
        probeset = lambda wc: f"--probe-set={config['probeset']}" if config.get("probeset") not in [False, "false", "False", None] else "",
        cytaimage  = lambda wc: f"--cyta-image={data_dir}/{wc.sample}_fastqs/{wc.sample}_cytaimage.jpeg" if config.get("cytaimage") not in [False, "false", "False", None] else "",
        loupealignment = lambda wc: f"--loupe-alignment={config['loupealignment']}" if config.get("loupealignment") not in [False, "false", None] else ""

    threads: 16

    resources:
        mem_mb=128000,
        runtime=3600,
        disc_mb=100000
    shell:
        """
        module load spaceranger
        rm -rf results/{wildcards.sample}

        spaceranger count \
            --id={wildcards.sample} \
            --output-dir=results/{wildcards.sample} \
            --transcriptome={params.transcriptome} \
            --fastqs={input.fastqs} \
            --image={input.image} \
            --slide={params.slide} \
            --area={params.area} \
            --create-bam={params.create_bam} \
            {params.nucleus_segmentation} \
            {params.loupealignment} \
            {params.probeset} \
            {params.cytaimage}
        """

rule seurat_process:
    input:
        matrix = get_seurat_input
    output:
        rds = "results/seurat/{sample}/seurat.rds",
        violin = "results/seurat/{sample}/qc_violin.pdf",
        umap = "results/seurat/{sample}/umap_plot.pdf",
        spatial = "results/seurat/{sample}/spatial_plot.pdf"
    resources:
        mem_mb=32000
    params:
        platform = lambda wc: PLATFORM_DICT[wc.sample]
    shell:
        """
        module load R
        mkdir -p results/seurat/{wildcards.sample}

        Rscript scripts/seurat_process.R \
        {wildcards.sample} \
        {params.platform} \
        {input.matrix} \
        {output.rds} \
        {output.violin} \
        {output.umap} \
        {output.spatial}
        """


#rule merge_seurat:
#    input:
#        samples_csv = "input/samples.csv",
#        seurat_objs = expand("results/{sample}/seurat.rds", sample=samples_df["sample"].tolist())
#    output:
#        rds = "results/merged_seurat.rds",
#        plot = "results/merged_umap.pdf",
#        de_condition = "results/de_condition.csv",
#        de_cluster = "results/de_by_cluster.csv"
#    shell:
#        """
#        Rscript scripts/merge_seurat.R {input.samples_csv} {output.rds} {output.plot} {output.de_condition} {output.de_cluster}
#        """

