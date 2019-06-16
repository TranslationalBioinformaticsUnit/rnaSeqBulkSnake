def get_trimmed(wildcards):
    if samples.loc[wildcards.sample, "fq2"]:
        # paired-end sample
        return expand(config["outputDir"] + "trimmed/{sample}.{group}.fastq.gz",
                      sample=wildcards.sample, group=[1, 2])
    # single end sample
    return config["outputDir"] + "trimmed/{sample}.fastq.gz"

rule make_star_index:
    input:
        seq = config["ref"]["sequence"],
        annotation = config["ref"]["annotation"]
    output:
        index = directory(config["outputDir"] + "star/reference_index")
    log:
        config["outputDir"] + "logs/star/reference_indexing.log"
    params:
        star_params =config["params"]["star_index"],
        overhang = config["ref"]["star_overhang"]
    threads: 12
    shell:"""
        mkdir -p {output} && \
        STAR --runThreadN {threads} --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.seq} \
            --sjdbGTFfile {input.annotation} \
            --sjdbOverhang {params.overhang} {params.star_params}
    """

rule align:
    input:
        sample=get_trimmed,
        index=config["outputDir"] + "star/reference_index"
    output:
        # see STAR manual for additional output files
        config["outputDir"] + "star/{sample}/Aligned.sortedByCoord.out.bam",
        config["outputDir"] + "star/{sample}/ReadsPerGene.out.tab"
    log:
        config["outputDir"] + "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        annotation  =config["ref"]["annotation"],
        star_params =config["params"]["star"],
        outPrefix   =config["outputDir"] + "star/{sample}/"
    threads:
        24
    shell:"""
        STAR --runThreadN {threads} --readFilesIn {input.sample} --genomeDir {input.index} --quantMode GeneCounts --sjdbGTFfile {params.annotation}\
        --outSAMtype BAM Unsorted  --readFilesCommand gunzip -c --outFileNamePrefix {params.outPrefix}  {params.star_params} > {log} 2>&1 && \
        samtools sort -@ {threads}  -o {params.outPrefix}/Aligned.sortedByCoord.out.bam {params.outPrefix}/Aligned.out.bam && \
        rm {params.outPrefix}/Aligned.out.bam
    """

rule qualimap_qc:
    input:  config["outputDir"] + "star/{sample}/Aligned.sortedByCoord.out.bam",
    output: config["outputDir"] + 'star/post_mapping_qualimap/{sample}/qualimapReport.html',
    params:
        outdir=config["outputDir"] + 'star/post_mapping_qualimap/{sample}',
        gtf=config["ref"]["annotation"],
    log: config["outputDir"] + "logs/star/qualimap/{sample}_qualimap.log"
    shell:"""
        qualimap rnaseq -bam {input} -gtf {params.gtf} --outdir {params.outdir} --java-mem-size=8G >    {log}

        """

rule run_picardmetrics:
    input: config["outputDir"] + "star/{sample}/Aligned.sortedByCoord.out.bam",
    output: config["outputDir"] + "star/bam_metrics/{sample}.metrics"
    shell:"""
        picard CollectInsertSizeMetrics I={input} H={output}.insertsize.pdf O={output} USE_JDK_INFLATER=true USE_JDK_DEFLATER=true

        """

rule create_insertsize_tsv:
    input:  config["outputDir"] + 'star/bam_metrics/{sample}.metrics'
    output: config["outputDir"] + 'star/bam_metrics/{sample}.insertsizes.tsv'
    shell:"""
        python ../scripts/collect_picard_metrics.py {input} {output}

        """