rule spades:
    input:
        hr1 = os.path.join(config["reads"], "host_removed", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["reads"], "host_removed", "{sample}_2_hr.fastq.gz"),
    threads: 8
    conda: "../envs/spades_env.yml"
    output:
        dir = directory(os.path.join(config["outdir"], "{sample}", "assembly")),
        contigs = os.path.join(config["outdir"], "{sample}", "assembly", "contigs.fasta")
    log:
        "logs/spades/{sample}.log"
    shell:
        """
        spades.py --meta -m 400 -t {threads} -1 {input.hr1} -2 {input.hr2} -o {output.dir} &>> {log}
        """
