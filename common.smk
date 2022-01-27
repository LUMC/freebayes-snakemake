containers = {
    "cutadapt": "docker://quay.io/biocontainers/cutadapt:3.5--py36hc5360cc_0",
    "debian": "docker://debian:latest",
}
default = {"setting1": "common.smk", "setting2": "common.smk", "setting3": "common.smk"}


def get_outfile():
    return "outputfile.txt"


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


def rg_per_sample():
    """Yield (sample, readgroup) for every sample"""
    for sample in pep.samples:
        sname = sample["sample_name"]
        fastq = pep.sample_table.loc[sname, "forward"]
        # If there is only a single fastq file
        if isinstance(fastq, str):
            yield sname, "rg1"
        # If there are multiple fastq files
        else:
            for i in range(len(fastq)):
                yield sname, f"rg{i+1}"
