import pandas as pd

shell.prefix("source config.sh; set -euo pipefail;")

configfile: "config.yaml"
samples = pd.read_table("samples.tsv", index_col=0)

samples['fq1'] = list(config["readDataDir"] + samples.index.astype(str) + "_1.fastq.gz")
samples['fq2'] = list(config["readDataDir"] + samples.index.astype(str) + "_2.fastq.gz")

rule all:
    input:
        expand(config["outputDir"] + 'qc/trimmed/{sample}.1_fastqc.html', sample=samples.index),
        expand(config["outputDir"] + 'qc/trimmed/{sample}.2_fastqc.html', sample=samples.index),
        expand(config["outputDir"] + 'qc/raw/{sample}_1_fastqc.html', sample=samples.index),
        expand(config["outputDir"] + 'qc/raw/{sample}_2_fastqc.html', sample=samples.index),
        expand(config["outputDir"] + "results/diffexp/{contrast}.diffexp.tsv",
               contrast=config["diffexp"]["contrasts"]),
        expand(config["outputDir"] + "results/diffexp/{contrast}.edgeR_diffexp.tsv",
               contrast=config["diffexp"]["contrasts"]),
        expand(config["outputDir"] + "results/diffexp/{contrast}.voom_diffexp.tsv",
               contrast=config["diffexp"]["contrasts"]),
        config["outputDir"] + "results/pca.svg",
        config["outputDir"] + 'multiqc_report/multiqc_report.html'


include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
