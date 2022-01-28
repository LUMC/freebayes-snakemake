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
        vcf=expand(
            "{sample}/{sample}.phased.vcf.gz", sample=pep.sample_table["sample_name"]
        ),
        bam=expand(
            "{sample}/{sample}.phased.bam", sample=pep.sample_table["sample_name"]
        ),


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
        rg="@RG\\tID:{sample}-{readgroup}\\tSM:{sample}",
    log:
        bwa="log/{sample}_{readgroup}.align.bwa.log",
        sam="log/{sample}_{readgroup}.align.samtools.log",
    container:
        containers["bwa-mem2"]
    shell:
        """
        bwa-mem2 mem \
            {input.reference} \
            -R '{params.rg}' \
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
        tbi="{sample}/{sample}.vcf.gz.tbi",
    log:
        "log/{sample}_call_variants.txt",
    container:
        containers["freebayes"]
    shell:
        """
        freebayes \
            --fasta-reference {input.reference} \
            --bam {input.bam} 2> {log} | bgzip > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule phase_variants:
    input:
        bam=rules.markdup.output.bam,
        reference=config["reference"],
        vcf=rules.call_variants.output.vcf,
    output:
        vcf="{sample}/{sample}.phased.vcf.gz",
        tbi="{sample}/{sample}.phased.vcf.gz.tbi",
    log:
        "log/{sample}_phase_variants.txt",
    container:
        containers["whatshap"]
    shell:
        """
        whatshap phase \
            --reference {input.reference} \
            {input.vcf} \
            {input.bam} 2> {log} | bgzip > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule phase_reads:
    input:
        bam=rules.markdup.output.bam,
        reference=config["reference"],
        vcf=rules.phase_variants.output.vcf,
    output:
        bam="{sample}/{sample}.phased.bam",
    log:
        "log/{sample}_phase_reads.txt",
    container:
        containers["whatshap"]
    shell:
        """
        whatshap haplotag \
            --reference {input.reference} \
            --output {output.bam} \
            {input.vcf} \
            {input.bam} 2> {log}
        """
