import pandas as pd
import os

"""========================================================================="""
"""                                 Parameters                              """
"""========================================================================="""

# File locations
data_dir = '/data/'
work_dir = os.getcwd()
meta_data = os.path.join(work_dir, 'input/samples.csv')

# Load sample metadata
samples_df = pd.read_csv(meta_data)
xenium_samples = samples_df[samples_df["platform"] == "xenium"]["sample"].tolist()
visium_samples = samples_df[samples_df["platform"] == "visium"]["sample"].tolist()

sample_platform_dict = samples_df.set_index("sample")["platform"].to_dict()

"""========================================================================="""
"""                                  Workflow                               """
"""========================================================================="""

# Final targets for both platforms
rule all:
    input:
        expand("results/{sample}/outs/filtered_feature_bc_matrix.h5", sample=visium_samples) +
        expand("results/{sample}/xenium_processed.tsv", sample=xenium_samples)

# SpaceRanger count for Visium
rule spaceranger_count:
    input:
        fastqs = lambda wildcards: os.path.join(data_dir, wildcards.sample, "fastq"),
        image = lambda wildcards: os.path.join(data_dir, wildcards.sample, "image.tif"),
        slide = lambda wildcards: os.path.join(data_dir, wildcards.sample, f"{wildcards.sample}.json")
    output:
        "results/{sample}/outs/filtered_feature_bc_matrix.h5"
    params:
        id = "{sample}",
        sample = "{sample}",
        transcriptome = config.transcriptome
    threads: 8
    shell:
        """
        module load spaceranger/4.0.1
        spaceranger count \
            --id={params.id} \
            --transcriptome={config.transcriptome} \
            --fastqs={input.fastqs} \
            --sample={params.sample} \
            --image={input.image} \
            --slide={input.slide} \
            --localcores={threads} \
            --probe-set={config.probeset}
        """

# Xenium CSV processing
rule xenium_process:
    input:
        matrix = lambda wildcards: os.path.join(data_dir, wildcards.sample, "cell_feature_matrix.csv")
    output:
        processed = "results/{sample}/xenium_processed.tsv"
    shell:
        """
        awk -F',' 'BEGIN{{OFS="\\t"}} NR==1{{print $0}} NR>1{{print $0}}' {input.matrix} > {output.processed}
        """

rule seurat_process:
    input:
        matrix = lambda wildcards: (
            f"results/{wildcards.sample}/outs/filtered_feature_bc_matrix.h5"
            if sample_platform_dict[wildcards.sample] == "visium" else
            f"results/{wildcards.sample}/xenium_processed.tsv"
        )
    output:
        rds = "results/{sample}/seurat.rds",
        violin = "results/{sample}/qc_violin.pdf",
        umap = "results/{sample}/umap_plot.pdf",
        spatial = "results/{sample}/spatial_plot.pdf"
    params:
        platform = lambda wildcards: sample_platform_dict[wildcards.sample]
    shell:
        """
        Rscript scripts/seurat_process.R {wildcards.sample} {params.platform} {input.matrix} \
        {output.rds} {output.violin} {output.umap} {output.spatial}
        """


rule merge_seurat:
    input:
        samples_csv = "input/samples.csv",
        seurat_objs = expand("results/{sample}/seurat.rds", sample=samples_df["sample"].tolist())
    output:
        rds = "results/merged_seurat.rds",
        plot = "results/merged_umap.pdf",
        de_condition = "results/de_condition.csv",
        de_cluster = "results/de_by_cluster.csv"
    shell:
        """
        Rscript scripts/merge_seurat.R {input.samples_csv} {output.rds} {output.plot} {output.de_condition} {output.de_cluster}
        """

