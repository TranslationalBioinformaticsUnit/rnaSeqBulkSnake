
def sample_id_to_str(wildcards=None):
    return ' '.join(samples.index)


rule count_matrix:
    input:
        expand(config["outputDir"] + "star/{sample}/ReadsPerGene.out.tab", sample=samples.index)
    output:
        config["outputDir"] + "counts/all.tsv"
    params:
        samples=sample_id_to_str(None)
    script:
        "../scripts/count-matrix.py"


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        counts  = config["outputDir"] + "counts/all.tsv",
        samples = "samples.tsv"
    output:
        config["outputDir"] + "deseq2/all.rds"
    conda:
        "../envs/deseq2.yaml"
    log:
        config["outputDir"] + "logs/deseq2/init.log"
    threads: get_deseq2_threads(None)
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        config["outputDir"] + "deseq2/all.rds"
    output:
        config["outputDir"] + "results/pca.svg"
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        config["outputDir"] + "logs/pca.log"
    script:
        "../scripts/plot-pca.R"


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]


rule deseq2:
    input:
        config["outputDir"] + "deseq2/all.rds"
    output:
        table   = config["outputDir"] + "results/diffexp/{contrast}.diffexp.tsv",
        ma_plot = config["outputDir"] + "results/diffexp/{contrast}.ma-plot.pdf",
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        config["outputDir"] + "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"


rule voom:
    input:
        counts  = config["outputDir"] + "counts/all.tsv",
        samples = "samples.tsv"
    output:
        table   = config["outputDir"] + "results/diffexp/{contrast}.voom_diffexp.tsv",
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        config["outputDir"] + "logs/edgeR/{contrast}.voom_diffexp.log"
    threads: 12
    script:
        "../scripts/edgeR.R"

rule edgeR:
    input:
        counts  = config["outputDir"] + "counts/all.tsv",
        samples = "samples.tsv"
    output:
        table   = config["outputDir"] + "results/diffexp/{contrast}.edgeR_diffexp.tsv",
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        config["outputDir"] + "logs/edgeR/{contrast}.edgeR_diffexp.log"
    threads: 12
    script:
        "../scripts/edgeR.R"


rule run_multiqc:
    input:
        expand(config["outputDir"] + 'qc/raw/{sample_name}_1_fastqc.html', sample_name=samples.index),
        expand(config["outputDir"] + 'qc/raw/{sample_name}_2_fastqc.html', sample_name=samples.index),
        expand(config["outputDir"] + 'star/bam_metrics/{sample}.metrics', sample=samples.index),
        expand(config["outputDir"] + 'star/post_mapping_qualimap/{sample}/qualimapReport.html', sample=samples.index),
        expand(config["outputDir"] + 'star/{sample}/Aligned.sortedByCoord.out.bam', sample=samples.index),
    output:
        config["outputDir"] + 'multiqc_report/multiqc_report.html'
    params:
        in_dir = config["outputDir"],
        out_dir = config["outputDir"] + "multiqc_report"
    log:
        config["outputDir"] + "logs/multiqc/multiqc_report.log"
    shell:
        'export LC_ALL=en_US.UTF-8 && multiqc -f --outdir {params.out_dir} {params.in_dir}'

