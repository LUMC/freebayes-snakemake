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
        stats=expand(
            "{sample}/{sample}.phased.tsv", sample=pep.sample_table["sample_name"]
        ),
        multiqc="multiqc_report.html",


rule cutadapt:
    input:
        fin=get_forward,
        rin=get_reverse,
    output:
        fout=temp("{sample}/{readgroup}_R1.fastq.gz"),
        rout=temp("{sample}/{readgroup}_R2.fastq.gz"),
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
        bam=temp("{sample}/{readgroup}.sorted.bam"),
        bai=temp("{sample}/{readgroup}.sorted.bam.bai"),
    params:
        compression_level=1,
        rg="@RG\\tID:{sample}-{readgroup}\\tSM:{sample}",
    log:
        bwa="log/{sample}_{readgroup}.align.bwa.log",
        sam="log/{sample}_{readgroup}.align.samtools.log",
    container:
        containers["bwa-mem2"]
    threads: 11
    shell:
        """
        # WARNING: This works fine for 1 thread (testing) and 11 threads, but
        # will under- or overprovision the tasks when other values are chosen

        # Use at most 8 threads for bwa-mem2
        bwa_cores=$(({threads} >= 8 ? 8 : 1))

        # Use at most 3 threads for samtools sort, but only when {threads} is
        # at least 8
        sam_cores=$(({threads} >= 8 ? 3 : 1))

        bwa-mem2 mem \
            {input.reference} \
            -R '{params.rg}' \
            -t $bwa_cores \
            {input.fin} {input.rin} 2> {log.bwa} |
            samtools sort \
                -l {params.compression_level} \
                -@ $sam_cores \
            - -o {output.bam} 2> {log.sam};
            samtools index {output}
        """


rule markdup:
    input:
        bam=get_bamfiles,
        bai=get_baifiles,
    output:
        bam=temp("{sample}/{sample}.bam"),
        bai=temp("{sample}/{sample}.bam.bai"),
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


rule collect_hs_metrics:
    input:
        bam=rules.markdup.output.bam,
        bai=rules.markdup.output.bai,
        bed=config.get("capture_bed", ""),
        reference=config["reference"],
        dictionary=config["reference"].split(".")[0] + ".dict",
    output:
        metrics="{sample}/hs_metrics.txt",
    log:
        "log/{sample}_hs_metrics.txt",
    container:
        containers["picard"]
    shell:
        """
        picard BedToIntervalList \
            INPUT={input.bed} \
            SEQUENCE_DICTIONARY={input.dictionary} \
            OUTPUT=capture.list

        picard CollectHsMetrics \
            I={input.bam} \
            R={input.reference} \
            BAIT_INTERVALS=capture.list \
            TARGET_INTERVALS=capture.list \
            O={output.metrics}
        """


rule call_variants:
    input:
        bam=rules.markdup.output.bam,
        bai=rules.markdup.output.bai,
        reference=config["reference"],
    output:
        vcf=temp("{sample}/{sample}.vcf.gz"),
        tbi=temp("{sample}/{sample}.vcf.gz.tbi"),
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
        tbi=rules.call_variants.output.tbi,
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
        bai=rules.markdup.output.bai,
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


rule phase_stats:
    input:
        vcf=rules.phase_variants.output.vcf,
    output:
        block="{sample}/{sample}.phased.blocklist",
        gtf="{sample}/{sample}.phased.gtf",
        tsv="{sample}/{sample}.phased.tsv",
    log:
        "log/{sample}_phase_stats.txt",
    container:
        containers["whatshap"]
    shell:
        """
        whatshap stats \
            --block-list {output.block} \
            --gtf {output.gtf} \
            --tsv {output.tsv} \
            {input.vcf}  2> {log}
        """


rule multiqc:
    input:
        markdup=expand(
            "log/{sample}_markdup.txt",
            sample=pep.sample_table["sample_name"],
        ),
        hs_metrics=get_hs_metrics,
        config=srcdir("config/multiqc_config.yml"),
    output:
        "multiqc_report.html",
    log:
        "log/multiqc.log",
    container:
        containers["multiqc"]
    shell:
        """
        rm -f multiqc_file_list.txt
        for file in {input.markdup} {input.hs_metrics}; do
            echo $file >> multiqc_file_list.txt
        done

        multiqc \
            --force \
            --file-list multiqc_file_list.txt \
            --config {input.config}
        """
