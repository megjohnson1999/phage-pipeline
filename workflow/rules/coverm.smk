rule rename_contigs:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}_DASTool_bins")
    output:
        dir = directory(os.path.join(config["outdir"], "all_bins")
        samples = directory(os.path.join(config["outdir"], "all_bins", "{sample}"))
    shell:
        """
        mkdir -p {output.dir}

        for file in {input}/*
        do
        sed "s/^>/>{wildcards.sample}_/" $file > {output.samples}_$(basename $file)
        done

        touch {output}
        """
        
rule coverm_cluster:
    input:
        os.path.join(config["outdir"], "all_bins")
    threads: 24
    conda: "../envs/coverm_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "all_bins_clustered"))
    shell:
        """
        coverm cluster -d {input} \
        --output-representative-fasta-directory-copy {output} \
        --threads {threads}

        # Concatenate clustered files
        cat {output}/* > {config[outdir]}/all_bins_clustered.fasta
        """

rule coverm_mapping:
    input:
        hr1 = os.path.join(config["reads"], "host_removed", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["reads"], "host_removed", "{sample}_2_hr.fastq.gz"),
        contigs = os.path.join(config["outdir"], "all_bins_clustered")
    threads: 24
    conda: "../envs/coverm_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "coverm"))
    log:
        os.path.join(config["outdir"], "logs", "coverm", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "coverm", "{sample}_bmrk.txt")
    shell:
        """
        mkdir -p {output}

        coverm contig -1 {input.hr1} -2 {input.hr2} -r {config[outdir]}/all_bins_clustered.fasta \
        --mapper minimap2-sr --threads {threads} \
        --methods rpkm count variance mean covered_fraction covered_bases \
        > {output}/{wildcards.sample}_stats.txt 2> {log}
        """
