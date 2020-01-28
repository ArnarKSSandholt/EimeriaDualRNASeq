# Snakemake pipeline for running a dual RNA-seq analysis on data from an Eimeria infection of chicken cells

configfile: "config.yaml"

rule trim:
    input:
        r1="data/{experiment}/{exp_meta}/{sample_group}/{sample}_R1_001.fastq.gz"
        r2="data/{experiment}/{exp_meta}/{sample_group}/{sample}_R2_001.fastq.gz"
    output:
        r1="results/trimmomatic/{experiment}/{sample}_R1_paired.fastq.gz"
        r2="results/trimmomatic/{experiment}/{sample}_R2_paired.fastq.gz"
        r1_unpaired="results/trimmomatic/{experiment}/{sample}_R1_unpaired.fastq.gz"
        r2_unpaired="results/trimmomatic/{experiment}/{sample}_R2_unpaired.fastq.gz"
    params:
        adapter=config["trim_params"]["adapter"]
        window=config["trim_params"]["window"]
    conda:
      "envs/trimmomatic.yaml"
    threads: 8
    shell:
        "trimmomatic PE -threads {threads} {input.r1} {input.r2}"
        "{output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired}"
        "{params.adapter} {params.window}"
    