from snakemake.utils import min_version
min_version("8.0")

configfile: "config.yaml"
SAMPLES, = glob_wildcards("reads/{sample}_1.fastq.gz")

include: "rules/preprocessing.smk"
include: "rules/assembly.smk"
include: "rules/binning.smk"
include: "rules/refine_bins.smk"
include: "rules/phages.smk"


rule all:
    input:
        expand("out/{sample}/phage_analysis/genomad", sample = SAMPLES)
    default_target: True
