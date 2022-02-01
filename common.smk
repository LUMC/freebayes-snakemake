containers = {
    # bwa-mem2 2.2.1, samtools 1.14
    "bwa-mem2": "docker://quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:bc6f1a049835e70d6fac927e979a8ad9819e1b56-0",
    "cutadapt": "docker://quay.io/biocontainers/cutadapt:3.5--py36hc5360cc_0",
    "debian": "docker://debian:latest",
    "freebayes": "docker://quay.io/biocontainers/freebayes:1.3.5--py36h74fc37f_4",
    "multiqc": "docker://quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0",
    "sambamba": "docker://quay.io/biocontainers/sambamba:0.8.1--hadffe2f_1",
    "whatshap": "docker://quay.io/biocontainers/whatshap:1.2.1--py38h8c62d01_0",
}

# Default settings for the pipeline
default = dict()


def get_forward(wildcards):
    return get_fastq(wildcards, "forward")


def get_reverse(wildcards):
    return get_fastq(wildcards, "reverse")


def get_fastq(wildcards, direction):
    fastq = pep.sample_table.loc[wildcards.sample, direction]

    # If a single fastq file is specified, we put it in a list
    if isinstance(fastq, str):
        fastq = [fastq]
        nr_fastq = 1
    # If multiple fastq files were specified, it is already a list
    else:
        nr_fastq = len(fastq)

    # Here, we use the readgroup wildcard to pick the correct fastq file
    readgroups = {f"rg{rg+1}": fastq for rg, fastq in zip(range(len(fastq)), fastq)}
    return readgroups[wildcards.readgroup]


def get_bamfiles(wildcards):
    """Return the bam files for a single sample"""
    sname = wildcards.sample
    fastq = pep.sample_table.loc[sname, "forward"]
    if isinstance(fastq, str):
        return f"{sname}/rg1.sorted.bam"
    else:
        return [f"{sname}/rg{i+1}.sorted.bam" for i in range(len(fastq))]


def get_baifiles(wildcards):
    """Return the bai files for a single sample"""
    bamfiles = get_bamfiles(wildcards)
    if isinstance(bamfiles, str):
        return f"{bamfiles}.bai"
    else:
        return [f"{bam}.bai" for bam in get_bamfiles(wildcards)]
