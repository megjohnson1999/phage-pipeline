# Trim adapters from raw reads
rule fastp:
    input:
        r1 = os.path.join(config["reads"],"{sample}_1.fastq.gz"),
        r2 = os.path.join(config["reads"],"{sample}_2.fastq.gz"),
    conda: "../envs/fastp_test.yaml"
    output:
        tr1 = os.path.join(config["reads"], "{sample}_1_trimmed.fastq.gz"),
        tr2 = os.path.join(config["reads"], "{sample}_2_trimmed.fastq.gz"),
    log:
        os.path.join(config["outdir"], "logs", "fastp", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "fastp", "{sample}_bmrk.txt")
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.tr1} -O {output.tr2} 2> {log}"

# Get database for host removal step
rule get_db:
    conda: "../envs/kneaddata.yaml"
    output:
        "ref/db_done"
    log:
        os.path.join(config["outdir"], "logs", "get_db")
    shell:
        """
        mkdir -p ref
        kneaddata_database --download human_genome bowtie2 ref 2> {log}
        touch ref/db_done
        """
    
# Remove host contamination
rule host_removal:
    input:
        tr1 = os.path.join(config["reads"], "{sample}_1_trimmed.fastq.gz"),
        tr2 = os.path.join(config["reads"], "{sample}_2_trimmed.fastq.gz"),
        db_done = "ref/db_done"
    params:
        db = config["host_database"]
    threads: 16
    conda: "../envs/bowtie2.yaml"
    output:
        hr1 = os.path.join(config["reads"], "host_removed", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["reads"], "host_removed", "{sample}_2_hr.fastq.gz"),
    log:
        os.path.join(config["outdir"], "logs", "host_removal", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "host_removal", "{sample}_bmrk.txt")
    shell:
        """
        mkdir -p {config[reads]}/host_removed
        bowtie2 -p {threads} -x {params.db} -1 {input.tr1} -2 {input.tr2} \
        --very-sensitive-local \
        --un-conc-gz {config[reads]}/host_removed/{wildcards.sample} \
        2> {log}
        mv {config[reads]}/host_removed/{wildcards.sample}.1 {output.hr1}
        mv {config[reads]}/host_removed/{wildcards.sample}.2 {output.hr2}
        """
