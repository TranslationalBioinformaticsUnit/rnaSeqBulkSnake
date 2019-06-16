#Bulk-RNAseq data processing for differential expression assessment

Gene-level differential expression assessment from bulk RNA-seq data with DESeq2, EdgeR, and Voom.

## How to use this pipeline

Step 1: Configure workflow
Set the input read file directory, the reference sequence, etc. in the config.yaml file.

Step 2: Provide sample annotation
Provide annotation of which sample corresponds to which treatment in the samples.tsv file.

Step 3: Setup your shell envinroment
Provide optional configuration of shell in config.sh (e.g. "module load <...>" or "export PATH=<...>" or "source activate <environment>")

Step 4: Test setup
Test configuration in dry-run: snakemake -n / snakemake --use-conda -n

Step 5:
Execute the workflow locally via

snakemake --use-conda --cores $ncpu
using $ncpu cores, or run it on a cluster via

runSnakemake.sh (with ressource requirement configuration in cluster.json)

## Workflow

- Read data QC with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Trimming with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- Alignment to reference with [STAR](https://github.com/alexdobin/STAR)
- Mapping QC with [qualimap](http://qualimap.bioinfo.cipf.es/) and [picard](https://broadinstitute.github.io/picard/)
- Differential expression assessment with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [Voom](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)
- QC report generation with [Multiqc](https://multiqc.info/)


## Test Dataset

- genome: ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.chromosome.21.fa.gz
- annotation: ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz
- reads: sampled from https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/HG00096.1.M_111124_6.bam

## Dependencies

To install:
- star
- multiqc
- fastqc
- trimmomatic
- qualimap
- picard

Via anaconda with --use-conda:
- deseq2
- biocparallel
- edger
- limma