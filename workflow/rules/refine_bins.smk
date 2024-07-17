rule dastool:
    input:
        contigs_filt = "out/{sample}/assembly/contigs_filt_1000bp.fasta",
        concoct_tsv = "out/{sample}/binning/dastool/concoct.contigs2bin.tsv",
        maxbin_tsv = "out/{sample}/binning/dastool/maxbin.contigs2bin.tsv",
        metabat_tsv = "out/{sample}/binning/dastool/metabat.contigs2bin.tsv"
    threads: config.get("num_threads", 8)
    conda: "../envs/dastool_env.yaml"
    output:
        "out/{sample}/binning/dastool/done"
    log:
        "logs/dastool/{sample}.log"
    shell:
        """
        touch {output}
        
        DAS_Tool --threads {threads} --write_bins --write_unbinned \
        -i {input.concoct_tsv},{input.maxbin_tsv},{input.metabat_tsv} \
        -l concoct,maxbin,metabat -c {input.contigs_filt} \
        -o out/{wildcards.sample}/binning/dastool/{wildcards.sample} \
        &>> {log} || true
        
        """

rule graphbin:
    input:
        contigs_filt = "out/{sample}/assembly/contigs_filt_1000bp.fasta",
        dastool = "out/{sample}/binning/dastool/done"
    threads: config.get("num_threads", 8)
    conda: "../envs/graphbin_env.yaml"
    output:
        directory("out/{sample}/binning/graphbin")
    log:
        "logs/graphbin/{sample}.log"
    shell:
        """
        # For samples with bins, format bins and run graphBin
        if [ -s out/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_contig2bin.tsv ]; then
        cat out/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_contig2bin.tsv \
        | tr -s "\\t" "," > out/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_contig2bin.csv

        python scripts/prepResult.py \
        --binned out/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins \
        --output "out/{wildcards.sample}/binning/dastool/" &>> {log}

        mkdir -p {output}

        graphbin --assembler spades \
        --graph out/{wildcards.sample}/assembly/assembly_graph_after_simplification.gfa \
        --contigs {input.contigs_filt} \
        --paths out/{wildcards.sample}/assembly/contigs.paths \
        --binned out/{wildcards.sample}/binning/dastool/initial_contig_bins.csv \
        --output {output}/{wildcards.sample} &>> {log}

        # For samples without bins, create necessary outputs
        else
        mkdir -p {output}
        fi

        """

rule checkm:
    input:
        "out/{sample}/binning/graphbin"
    threads: config.get("num_threads", 8)
    conda: "envs/checkm_env.yaml"
    output:
        directory("out/{sample}/binning/checkm")
    log:
        "logs/checkm/{sample}.log"
    shell:
        """
        mkdir -p {output}
        if [ -s {input}/{wildcards.sample}bins ]; then
        checkm lineage_wf -x fasta {input}/{wildcards.sample}bins/ \
        {output}/ -t {threads} --tab_table -f {output}/checkm_out.tsv &>> {log}
        fi
        """ 
