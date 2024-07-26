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
        os.path.join(config["logs"], "dastool", "{sample}.log")
    shell:
        """
        touch {output}
        
        DAS_Tool --threads {threads} --write_bins --write_unbinned \
        -i {input.concoct_tsv},{input.maxbin_tsv},{input.metabat_tsv} \
        -l concoct,maxbin,metabat -c {input.contigs_filt} \
        -o {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample} \
        &> {log} || true
        
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
        os.path.join(config["logs"], "graphbin", "{sample}.log")
    shell:
        """
        # For samples with bins, format bins and run graphBin
        if [ -s {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_contig2bin.tsv ]; then
        cat {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_contig2bin.tsv \
        | tr -s "\\t" "," > {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_contig2bin.csv

        python scripts/prepResult.py \
        --binned {config[outdir]}/{wildcards.sample}/binning/dastool/{wildcards.sample}_DASTool_bins \
        --output {config[outdir]}/{wildcards.sample}/binning/dastool/ &> {log}

        mkdir -p {output}

        graphbin --assembler spades \
        --graph {config[outdir]}/{wildcards.sample}/assembly/assembly_graph_after_simplification.gfa \
        --contigs {input.contigs_filt} \
        --paths {config[outdir]}/{wildcards.sample}/assembly/contigs.paths \
        --binned {config[outdir]}/{wildcards.sample}/binning/dastool/initial_contig_bins.csv \
        --output {output}/{wildcards.sample} &>> {log}

        # For samples without bins, create necessary outputs
        else
        mkdir -p {output}
        fi

        """

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

        # Extract the sequences >=4000bp, separate into individual fasta files, move to directory
        cat {config[outdir]}/{wildcards.sample}/binning/unbinned_contigs.fasta \
        | seqkit seq -m 4000 \
        > {config[outdir]}/{wildcards.sample}/binning/filt_4000_seqs_to_keep.fasta

        cd {config[outdir]}/{wildcards.sample}/binning
        cat filt_4000_seqs_to_keep.fasta \
        | awk '{\
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")}
        print $0 >> filename
        close(filename)
        }'
        mv NODE* graphbin/{wildcards.sample}bins
        cd ../../../..

        # Extract the IDs of the unbinned sequences <4000bp
        cat {config[outdir]}/{wildcards.sample}/binning/unbinned_contigs.fasta \
        | seqkit seq -n -M 3999 \
        > {config[outdir]}/{wildcards.sample}/binning/filt_4000_seqs_to_discard.txt
        # From main contigs file, get all sequences except for those on this list
        seqkit grep -v -f {config[outdir]}/{wildcards.sample}/binning/filt_4000_seqs_to_discard.txt \
        {input.contigs_filt} -o {output.final_contigs}

        mv {input.graphbin}/{wildcards.sample}bins/bin_unbinned.fa.fasta {input.graphbin}


        # Otherwise, if no bins were determined for the sample...
        else
        echo "No sequences were binned for sample {wildcards.sample}"
        # Filter sequences <4000bp from main contigs file
        cat {input.contigs_filt} | seqkit seq -m 4000 > {output.final_contigs}

        mkdir -p {input}/{wildcards.sample}bins
        cd {config[outdir]}/{wildcards.sample}/binning
        cat final_filtered_contigs.fasta \
        | awk '{\
        if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")}
        print $0 >> filename
        close(filename)
        }'
        mv NODE* graphbin/{wildcards.sample}bins
        cd ../../../..

        fi

        # Filter contigs for phispy input (5000bp filter)
        cat {output.final_contigs} | seqkit seq -m 5000 > {output.contigs_5000bp}
        """

rule checkm:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "graphbin")
    threads: 24
    conda: "../envs/checkm_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "binning", "checkm"))
    log:
        os.path.join(config["logs"], "checkm", "{sample}.log")
    shell:
        """
        mkdir -p {output}
        checkm lineage_wf -x fasta {input}/{wildcards.sample}bins/ \
        {output}/ -t {threads} --tab_table -f {output}/checkm_out.tsv &> {log}
        """ 
