rule filter_unbinned:
    input:
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta"),
        graphbin = os.path.join(config["outdir"], "{sample}", "binning", "graphbin")
    conda: "../envs/bowtie2.yaml"
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
    output: directory("ref/genomad_db")
    conda: "../envs/genomad_env.yaml"
    shell:
        "genomad download-database ref"

rule genomad:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        db = "ref/genomad_db"
    threads: config.get("num_threads", 8)
    conda: "../envs/genomad_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"))
    log:
        "logs/genomad/{sample}.log"
    shell:
        """
        # Create the output directory
        mkdir -p {output}
        # Run genomad
        genomad end-to-end --cleanup --threads {threads} \
        {input.contigs} {output} {input.db} &> {log}
        """

rule download_bakta_db:
    output: directory("ref/bakta_db")
    conda: "../envs/bakta_env.yaml"
    log:
        "logs/download_bakta_db"
    shell:
        """
        bakta_db download --output {output} --type full \
        &> {log}
        """   

rule bakta:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filt_contigs_5000.fasta"),
        db = "ref/bakta_db"
    threads: config.get("num_threads", 8)
    conda: "../envs/bakta_env.yaml"
    output: 
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "bakta"))
    log:
        "logs/bakta/{sample}.log"
    shell:
        """
        bakta --db {input.db}/db --force --output {output} \
        --threads {threads} {input.contigs} &> {log}
        """

rule phispy:
    input:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "bakta")
    threads: config.get("num_threads", 8)
    conda: "../envs/phispy_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy"))
    log:
        "logs/phispy/{sample}.log"
    shell:
        "PhiSpy.py {input}/*.gbff -o {output} --output_choice 63 &> {log} || true"

rule phage_all:
    input:
        genomad = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta")
        phispy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy")
    output:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "done")
    shell:
        "echo 'Done'"
