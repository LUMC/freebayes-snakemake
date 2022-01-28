include: "common.smk"


pepfile: config["pepfile"]


# Apply the settings from the pepfile, overwriting the default ones
default.update(pep.config.get("freebayes-snakemake", dict()))

# Apply the options specified to snakemake, overwriting the default settings
# and the settings from the PEP file
default.update(config)

# Set the updated dict as the configuration for the pipeline
config = default


rule all:
    input:
        outfile=get_outfile(),
        samples=expand("{sample}.txt", sample=pep.sample_table["sample_name"]),
        trimmed=[f"{sample}/{rg}_R1.fastq.gz" for sample, rg in rg_per_sample()],
        bams=expand("{sample}/{sample}.bam", sample=pep.sample_table["sample_name"]),
        vcf=expand("{sample}/{sample}.vcf.gz", sample=pep.sample_table["sample_name"]),
        settings="settings.txt",


rule cutadapt:
    input:
        fin=get_forward,
        rin=get_reverse,
    output:
        fout="{sample}/{readgroup}_R1.fastq.gz",
        rout="{sample}/{readgroup}_R2.fastq.gz",
    log:
        "log/{sample}_{readgroup}_cutadapt.txt",
    container:
        containers["cutadapt"]
    threads: 4
    shell:
        """
        cutadapt \
            -a AGATCGGAAGAG \
            -A AGATCGGAAGAG \
            --compression-level=1 \
            --cores {threads} \
            --output {output.fout} \
            --paired-output {output.rout} \
            {input.fin} {input.rin} \
            > {log}
        """


rule align:
    input:
        fin=rules.cutadapt.output.fout,
        rin=rules.cutadapt.output.rout,
        reference=config["reference"],
    output:
        bam="{sample}/{readgroup}.sorted.bam",
        bai="{sample}/{readgroup}.sorted.bam.bai",
    params:
        compression_level=1,
    log:
        bwa="log/{sample}_{readgroup}.align.bwa.log",
        sam="log/{sample}_{readgroup}.align.samtools.log",
    container:
        containers["bwa-mem2"]
    shell:
        """
        bwa-mem2 mem \
            {input.reference} \
            {input.fin} {input.rin} 2> {log.bwa} |
            samtools sort -l {params.compression_level} \
            - -o {output.bam} 2> {log.sam};
            samtools index {output}
        """


rule markdup:
    input:
        bam=get_bamfiles,
    output:
        bam="{sample}/{sample}.bam",
        bai="{sample}/{sample}.bam.bai",
    params:
        compression_level=1,
    log:
        "log/{sample}_markdup.txt",
    container:
        containers["sambamba"]
    threads: 4
    shell:
        """
        sambamba markdup \
            --nthreads={threads} \
            --compression-level={params.compression_level} \
            {input.bam} {output.bam} 2> {log}
        """


rule call_variants:
    input:
        bam=rules.markdup.output.bam,
        reference=config["reference"],
    output:
        vcf="{sample}/{sample}.vcf.gz",
    log:
        "log/{sample}_call_variants.txt",
    container:
        containers["freebayes"]
    shell:
        """
        freebayes \
            --fasta-reference {input.reference} \
            --bam {input.bam} | bgzip > {output.vcf} 2> {log}
        """


rule example:
    output:
        get_outfile(),
    log:
        "log/stdout.txt",
    container:
        containers["debian"]
    shell:
        """
        echo "Hello world!" > {output} 2> {log}
        """


rule sample:
    output:
        "{sample}.txt",
    log:
        "log/{sample}_touch.txt",
    container:
        containers["debian"]
    shell:
        """
        touch {output} 2> {log}
        """


rule settings:
    output:
        "settings.txt",
    params:
        s1=config["setting1"],
        s2=config["setting2"],
        s3=config["setting3"],
    log:
        "log/settings.txt",
    container:
        containers["debian"]
    shell:
        """
        echo {params.s1} {params.s2} {params.s3} > {output}
        """
