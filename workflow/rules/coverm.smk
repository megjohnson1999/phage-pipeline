rule rename_contigs:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta")
    output:
        os.path.join(config["outdir"], "{sample}", "binning", "renamed_final_filtered_contigs.fasta")
    shell:
        "sed "s/^>/>{wildcards.sample}_/" {input} > {output}"

rule concatenate_contigs:
    input:
        expand()
    output:
        os.path.join(config["outdir"], "concatenated_contigs.fasta")
    shell:
        """

rule cluster_and_derep:
    input:
        os.path.join(config["outdir"], "concatenated_contigs.fasta")
    output:
        os.path.join(config["outdir"], "dereplicated_contigs.fasta")
    shell:

rule coverm:
