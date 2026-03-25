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
        probeset = lambda wc: f"--probe-set={config['probeset']}" if config.get("probeset") not in [False, "false", "False", None] else ""
        cytaimage  = lambda wc: f"--cyta-image={data_dir}/{wc.sample}_fastqs/{wc.sample}_cytaimage.jpeg" if config.get("cytaimage") not in [False, "false", "False", None] else ""

    threads: 16

    resources:
        mem_mb=128000,
        runtime=3600,
        disc_mb=100000
    shell:
        """
        rm -rf results/{wildcards.sample}

        spaceranger count \
            --id={wildcards.sample} \
            --output-dir=results/{wildcards.sample} \
            --transcriptome={params.transcriptome} \
            --fastqs={input.fastqs} \
            --image={input.image} \
            --cytaimage={params.cytaimage} \
            --slide={params.slide} \
            --area={params.area} \
            --create-bam={params.create_bam} \
            --nucleus-segmentation={params.nucleus_segmentation} \
            {params.probeset}
        """

# Xenium CSV processing
#rule xenium_process:
#    input:
#        matrix = lambda wildcards: os.path.join(data_dir, wildcards.sample, "cell_feature_matrix.csv")
#    output:
#        processed = "results/{sample}/xenium_processed.tsv"
#    shell:
#        """
#        awk -F',' 'BEGIN{{OFS="\\t"}} NR==1{{print $0}} NR>1{{print $0}}' {input.matrix} > {output.processed}
#        """

rule seurat_process:
    input:
        matrix = "results/{sample}/outs/"
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

