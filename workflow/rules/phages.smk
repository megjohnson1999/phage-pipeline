rule genomad_db:
    output: directory("ref/genomad_db")
    conda: "../envs/genomad_env.yaml"
    shell:
        "genomad download-database ref"

rule genomad:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta"),
        db = "ref/genomad_db"
    threads: 24
    conda: "../envs/genomad_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"))
    log:
        os.path.join(config["logs"], "genomad", "{sample}.log")
    shell:
        """
        # Create the output directory
        mkdir -p {output}
        # Run genomad
        genomad end-to-end --cleanup --threads {threads} \
        --splits 16 \
        {input.contigs} {output} {input.db} &> {log}
        """

rule download_bakta_db:
    output: directory("ref/bakta_db")
    conda: "../envs/bakta_env.yaml"
    log:
        os.path.join(config["logs"], "download_bakta_db")
    shell:
        """
        bakta_db download --output {output} --type full \
        &> {log}
        """   

rule bakta:
    input:
        contigs = os.path.join(config["outdir"], "{sample}", "binning", "final_filt_contigs_5000.fasta"),
        db = "ref/bakta_db"
    threads: 24
    conda: "../envs/bakta_env.yaml"
    output: 
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "bakta"))
    log:
        os.path.join(config["logs"], "bakta", "{sample}.log")
    shell:
        """
        bakta --db {input.db}/db --force --output {output} \
        --threads {threads} {input.contigs} &> {log}
        """

rule phispy:
    input:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "bakta")
    threads: 24
    conda: "../envs/phispy_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy"))
    log:
        os.path.join(config["logs"], "phispy", "{sample}.log")
    shell:
        "PhiSpy.py {input}/*.gbff -o {output} --output_choice 63 &> {log} || true"

rule phage_all:
    input:
        checkm = os.path.join(config["outdir"], "{sample}", "binning", "checkm"),
        genomad = os.path.join(config["outdir"], "{sample}", "phage_analysis", "genomad"),
        phispy = os.path.join(config["outdir"], "{sample}", "phage_analysis", "phispy")
    output:
        os.path.join(config["outdir"], "{sample}", "phage_analysis", "done")
    shell:
        "touch {config[outdir]}/{wildcards.sample}/phage_analysis/done"
