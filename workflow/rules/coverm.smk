rule rename_contigs:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta")
    output:
        os.path.join(config["outdir"], "{sample}", "binning", "renamed_final_filtered_contigs.fasta")
    shell:
        'sed "s/^>/>{wildcards.sample}_/" {input} > {output}'

rule concatenate_contigs:
    input:
        expand(os.path.join(config["outdir"], "{sample}", "binning", "renamed_final_filtered_contigs.fasta"), sample = SAMPLES)
    output:
        os.path.join(config["outdir"], "concatenated_contigs.fasta")
    shell:
        'cat {input} > {output}'

rule cluster_and_derep:
    input:
        os.path.join(config["outdir"], "concatenated_contigs.fasta")
    threads: 24
    conda: "../envs/cdhit_env.yaml"
    output:
        os.path.join(config["outdir"], "dereplicated_contigs.fasta")
    shell:
        "cd-hit-est -i {input} -o {output} -c 0.95 -n 10 -T {threads}"

rule coverm_mapping:
    input:
        hr1 = os.path.join(config["reads"], "host_removed", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["reads"], "host_removed", "{sample}_2_hr.fastq.gz"),
        contigs = os.path.join(config["outdir"], "dereplicated_contigs.fasta")
    threads: 24
    conda: "../envs/coverm_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "coverm"))
    shell:
        """
        mkdir -p {output}

        coverm contig -1 {input.hr1} -2 {input.hr2} -r {input.contigs} \
        --mapper minimap2-sr --threads {threads} \
        --methods rpkm count variance mean covered_fraction covered_bases \
        > {output}/{wildcards.sample}_stats.txt
        """
