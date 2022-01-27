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
        bams=[f"{sample}_{rg}.bam" for sample, rg in rg_per_sample()],
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


rule map:
    input:
        f=get_forward,
        r=get_reverse,
    output:
        "{sample}_{readgroup}.bam",
    log:
        "log/{sample}_{readgroup}map.txt",
    container:
        containers["debian"]
    shell:
        """
        echo mem ref.fa {input.f} {input.r} > {output}
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
