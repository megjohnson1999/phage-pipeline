rule filter_unbinned:
    input:
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta"),
        graphbin = os.path.join(config["outdir"], "{sample}", "binning", "graphbin")
    conda: "../envs/minimap_env.yaml"
    output:
        final_contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        contigs_5000bp = os.path.join(config["outdir"], "{sample}", "binning", "final_filt_contigs_5000.fasta")
    shell:
        """
        # If there were bins for the sample...
        if [ -s {input.graphbin}/{wildcards.sample}graphbin_unbinned.csv ]; then
        # Get fasta file of unbinned sequences
        seqkit grep -f {input.graphbin}/{wildcards.sample}graphbin_unbinned.csv \
        {input.contigs_filt} -o {config[outdir]}/{wildcards.sample}/binning/unbinned_contigs.fasta
        # Extract the IDs of the unbinned sequences <4000bp
        cat {config[outdir]}/{wildcards.sample}/binning/unbinned_contigs.fasta \
        | seqkit seq -n -M 4000 \
        > {config[outdir]}/{wildcards.sample}/binning/filt_4000_seqs_to_discard.txt
        # From main contigs file, get all sequences except for those on this list
        seqkit grep -v -f {config[outdir]}/{wildcards.sample}/binning/filt_4000_seqs_to_discard.txt \
        {input.contigs_filt} -o {output.final_contigs}

        # Otherwise, if no bins were determined for the sample...
        else
        echo "No sequences were binned for sample {wildcards.sample}"
        # Filter sequences <4000bp from main contigs file
        cat {input.contigs_filt} | seqkit seq -m 4000 > {output.final_contigs}
        fi

        # Filter contigs for phispy input (5000bp filter)
        cat {output.final_contigs} | seqkit seq -m 5000 > {output.contigs_5000bp}
        """

rule genomad_db:
    output: directory(config["genomad_database"])
    conda: "../envs/genomad_env.yaml"
    shell:
        "genomad download-database ref"

rule genomad:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        db = config["genomad_database"]
    threads: 24
    conda: "../envs/genomad_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"))
    log:
        os.path.join(config["outdir"], "logs", "genomad", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "genomad", "{sample}_bmrk.txt")
    shell:
        """
        # Create the output directory
        mkdir -p {output}
        # Run genomad
        genomad end-to-end --cleanup --threads {threads} \
        --splits 16 \
        {input.contigs} {output} {input.db} 2> {log}
        """

rule download_bakta_db:
    output: directory(config["bakta_database"])
    conda: "../envs/bakta_env.yaml"
    shell:
        "bakta_db download --output {output} --type full"

rule bakta:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filt_contigs_5000.fasta"),
        db = config["bakta_database"]
    threads: 24
    conda: "../envs/bakta_env.yaml"
    output: 
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "bakta"))
    log:
        os.path.join(config["outdir"], "logs", "bakta", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "bakta", "{sample}_bmrk.txt")
    shell:
        """
        bakta --db {input.db}/db --force --skip-plot --output {output} \
        --threads {threads} {input.contigs} 2> {log}
        """

rule phispy:
    input:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "bakta")
    threads: 24
    conda: "../envs/phispy_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy"))
    log:
        os.path.join(config["outdir"], "logs", "phispy", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "phispy", "{sample}_bmrk.txt")
    shell:
        "PhiSpy.py {input}/*.gbff -o {output} --output_choice 63 2> {log} || true"

rule phage_all:
    input:
        cat = os.path.join(config["outdir"], "{sample}", "taxonomy", "CAT"),
        checkm = os.path.join(config["outdir"], "{sample}", "binning", "checkm"),
        genomad = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"),
        phispy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy")
    output:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "done")
    shell:
        "touch {config[outdir]}/{wildcards.sample}/phage_analysis/done"
