def get_fastq(wildcards):
    return samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()


rule raw_qc:
    input:
        r1=config["readDataDir"]+"{sample}_1.fastq.gz",
        r2=config["readDataDir"]+"{sample}_2.fastq.gz"
    params:
        out_dir = config["outputDir"] + 'qc/raw'
    output:
       config["outputDir"] + 'qc/raw/{sample}_1_fastqc.html',
       config["outputDir"] + 'qc/raw/{sample}_1_fastqc.zip',
       config["outputDir"] + 'qc/raw/{sample}_2_fastqc.html',
       config["outputDir"] + 'qc/raw/{sample}_2_fastqc.zip',
    shell:"""
            fastqc -o {params.out_dir} -f fastq {input.r1} {input.r2}
        """

rule trimmed_qc:
    input:
        r1=config["outputDir"] + "trimmed/{sample}.1.fastq.gz",
        r2=config["outputDir"] + "trimmed/{sample}.2.fastq.gz"
    params:
        out_dir = config["outputDir"] + 'qc/trimmed'
    output:
       config["outputDir"] + 'qc/trimmed/{sample}.1_fastqc.html',
       config["outputDir"] + 'qc/trimmed/{sample}.1_fastqc.zip',
       config["outputDir"] + 'qc/trimmed/{sample}.2_fastqc.html',
       config["outputDir"] + 'qc/trimmed/{sample}.2_fastqc.zip',
    shell:"""
            fastqc -o {params.out_dir} -f fastq {input.r1} {input.r2}
        """


rule trimmomatic_pe:
    input:
        r1 = config["readDataDir"] + "{sample}_1.fastq.gz",
        r2 = config["readDataDir"] + "{sample}_2.fastq.gz"
    output:
        r1 = config["outputDir"] + "trimmed/{sample}.1.fastq.gz",
        r2 = config["outputDir"] + "trimmed/{sample}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired = config["outputDir"] + "trimmed/{sample}.1.unpaired.fastq.gz",
        r2_unpaired = config["outputDir"] + "trimmed/{sample}.2.unpaired.fastq.gz"
    log:
        config["outputDir"] + "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:" + config["params"]["trimmo_adapter_path"] + ":2:30:10", "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"],
        # optional parameters
        extra="-phred33 -threads " + config["params"]["trimmo_threads"]
    wrapper:
        "0.17.4/bio/trimmomatic/pe"
