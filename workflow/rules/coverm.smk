rule rename_contigs:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "dastool", "{sample}.bins")
    output:
        directory(os.path.join(config["outdir"], "all_bins", "{sample}"))
    shell:
        """
        mkdir -p {output}
        files={config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins/*

        for file in $files
        do
        echo $file
        sed "s/^>/>{wildcards.sample}_/" $file \
        > {output}/{wildcards.sample}_$(basename $file)
        done
        """

rule move_files:
    input:
        expand(os.path.join(config["outdir"], "all_bins", "{sample}"), sample=SAMPLES)
    output:
        directory(os.path.join(config["outdir"], "all_bins", "all_samples"))
    shell:
        """
        mkdir -p {output}

        for dir in {input}
        do
        if [ -d "$dir" ]; then
        mv "$dir"/* {output}
        fi
        done
        """

        
rule coverm_cluster:
    input:
        os.path.join(config["outdir"], "all_bins", "all_samples")
    threads: 24
    conda: "../envs/coverm_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "all_bins_clustered"))
    log:
        os.path.join(config["outdir"], "logs", "coverm_cluster.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "coverm_cluster_bmrk.txt")
    shell:
        """
        coverm cluster --genome-fasta-directory {input} \
        -x .fa \
        --output-representative-fasta-directory-copy {output} \
        --threads {threads}

        # Concatenate clustered files
        cat {output}/* > {config[outdir]}/all_bins_clustered.fasta
        """

rule coverm_mapping:
    input:
        hr1 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_2_hr.fastq.gz"),
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
