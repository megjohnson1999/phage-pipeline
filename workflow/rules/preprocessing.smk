# Trim adapters from raw reads
rule fastp:
    input:
        r1 = os.path.join(config["reads"], config["fastq_names_1"]),
        r2 = os.path.join(config["reads"], config["fastq_names_2"]),
    params:
        l = config["fastp_min_sequence_length"]
    conda: "../envs/fastp_test.yaml"
    output:
        tr1 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_1_trimmed.fastq.gz"),
        tr2 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_2_trimmed.fastq.gz"),
    log:
        os.path.join(config["outdir"], "logs", "fastp", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "fastp", "{sample}_bmrk.txt")
    shell:
        "fastp -l {params.l} -i {input.r1} -I {input.r2} -o {output.tr1} -O {output.tr2} 2> {log}"

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
        tr1 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_1_trimmed.fastq.gz"),
        tr2 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_2_trimmed.fastq.gz"),
        db_done = "ref/db_done"
    params:
        db = config["human_ref"]
    threads: 16
    conda: "../envs/minimap_env.yaml"
    output:
        hr1 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_2_hr.fastq.gz"),
    log:
        os.path.join(config["outdir"], "logs", "host_removal", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "host_removal", "{sample}_bmrk.txt")
    shell:
        """
        minimap2 -ax sr {params.db} {input.tr1} {input.tr2} \
        | samtools view -bh \
        | samtools sort -o {config[outdir]}/{wildcards.sample}/preprocessing/{wildcards.sample}_output.bam
        samtools index {config[outdir]}/{wildcards.sample}/preprocessing/{wildcards.sample}_output.bam
        # Use samtools to get the reads that didn't map to host
        samtools fastq -F 3584 -f 77 {config[outdir]}/{wildcards.sample}/preprocessing/{wildcards.sample}_output.bam  \
        | gzip -c > {output.hr1}
        samtools fastq -F 3584 -f 141 {config[outdir]}/{wildcards.sample}/preprocessing/{wildcards.sample}_output.bam \
        | gzip -c > {output.hr2}
        """
