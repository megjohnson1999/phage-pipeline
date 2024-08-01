rule dastool:
    input: 
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta"),
        concoct_tsv = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "concoct.contigs2bin.tsv"),
        maxbin_tsv = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "maxbin.contigs2bin.tsv"),
        metabat_tsv = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "metabat.contigs2bin.tsv")
    threads: 24
    conda: "../envs/dastool_env.yaml"
    output:
        os.path.join(config["outdir"], "{sample}", "binning", "done")
    log:
        os.path.join(config["outdir"], "logs", "dastool", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "dastool", "{sample}_bmrk.txt")
    shell:
        """
        touch {output}
        
        DAS_Tool --threads {threads} --write_bins --write_unbinned \
        -i {input.concoct_tsv},{input.maxbin_tsv},{input.metabat_tsv} \
        -l concoct,maxbin,metabat -c {input.contigs_filt} \
        -o {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample} \
        2> {log} || true
        
        """

rule graphbin:
    input:
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta"),
        dastool = os.path.join(config["outdir"], "{sample}", "binning", "done")
    threads: 24
    conda: "../envs/graphbin_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "binning", "graphbin"))
    log:
        os.path.join(config["outdir"], "logs", "graphbin", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "graphbin", "{sample}_bmrk.txt")
    shell:
        """
        # For samples with bins, format bins and run graphBin
        if [ -s {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_contig2bin.tsv ]; then
        cat {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_contig2bin.tsv \
        | tr -s "\\t" "," > {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_contig2bin.csv

        python scripts/prepResult.py \
        --binned {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins \
        --output {config[outdir]}/{wildcards.sample}/binning/dastool/ 2> {log}

        mkdir -p {output}

        graphbin --assembler spades \
        --graph {config[outdir]}/{wildcards.sample}/assembly/assembly_graph_after_simplification.gfa \
        --contigs {input.contigs_filt} \
        --paths {config[outdir]}/{wildcards.sample}/assembly/contigs.paths \
        --binned {config[outdir]}/{wildcards.sample}/binning/dastool/initial_contig_bins.csv \
        --output {output}/{wildcards.sample} 2>> {log}

        # For samples without bins, create necessary outputs
        else
        mkdir -p {output}
        fi

        """

rule checkm:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "graphbin")
    threads: 24
    conda: "../envs/checkm_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "binning", "checkm"))
    log:
        os.path.join(config["outdir"], "logs", "checkm", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "checkm", "{sample}_bmrk.txt")
    shell:
        """
        mkdir -p {output}
        if [ -s {input}/{wildcards.sample}bins ]; then
        checkm lineage_wf -x fasta {input}/{wildcards.sample}bins/ \
        {output}/ -t {threads} --tab_table -f {output}/checkm_out.tsv 2> {log}
        fi
        """ 
