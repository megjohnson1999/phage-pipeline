rule binning_prep:
    input:
        hr1 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_2_hr.fastq.gz"),
        contigs = os.path.join(config["outdir"], "{sample}", "assembly", "contigs.fasta")
    threads: 12
    conda: "../envs/minimap_env.yaml"
    output:
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta"),
        sam = os.path.join(config["outdir"], "{sample}", "binning", "map_reads", "{sample}.sam"),
        bam = os.path.join(config["outdir"], "{sample}", "binning", "map_reads", "{sample}.bam"),
        sorted_bam = os.path.join(config["outdir"], "{sample}", "binning", "map_reads", "sorted_bam", "{sample}_sorted.bam")
    log:
        os.path.join(config["outdir"], "logs", "binning_prep", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "binning_prep", "{sample}_bmrk.txt")
    shell:
        """
        cat {input.contigs} | seqkit seq -m 1000 > {output.contigs_filt}

        bowtie2-build {output.contigs_filt} \
        {config[outdir]}/{wildcards.sample}/binning/map_reads/btdb 2> {log}

        bowtie2 --no-unal -p {threads} -x {config[outdir]}/{wildcards.sample}/binning/map_reads/btdb \
        -1 {input.hr1} -2 {input.hr2} -S {output.sam} 2>> {log}

        samtools view -@ {threads} -Sb -o {output.bam} {output.sam} 2>> {log}

        samtools sort -O bam -o {output.sorted_bam} {output.bam} 2>> {log}

        samtools index {output.sorted_bam} 2>> {log}
        """

rule concoct:
    input:
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta"),
        sorted_bam = os.path.join(config["outdir"], "{sample}", "binning", "map_reads", "sorted_bam", "{sample}_sorted.bam")
    threads: 24
    conda: "../envs/concoct_env.yaml"
    output:
        bed = os.path.join(config["outdir"], "{sample}", "binning", "concoct.out", "contigs_10k.bed"),
        fa = os.path.join(config["outdir"], "{sample}", "binning", "concoct.out", "contigs_10k.fasta"),
        cov_table = os.path.join(config["outdir"], "{sample}", "binning", "concoct.out", "coverage_table.tsv"),
        merged_csv = os.path.join(config["outdir"], "{sample}", "binning", "concoct.out", "results", "clustering_merged.csv")
    log:
        os.path.join(config["outdir"], "logs", "concoct", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "concoct", "{sample}_bmrk.txt")
    shell:
        """
        mkdir -p {config[outdir]}/{wildcards.sample}/binning/concoct.out
        mkdir -p {config[outdir]}/{wildcards.sample}/binning/concoct.out/results

        cut_up_fasta.py {input.contigs_filt} -c 10000 -o 0 --merge_last -b {output.bed} > {output.fa}

        concoct_coverage_table.py \
        {config[outdir]}/{wildcards.sample}/binning/concoct.out/contigs_10k.bed \
        {input.sorted_bam} > {output.cov_table}

        concoct --composition_file {config[outdir]}/{wildcards.sample}/binning/concoct.out/contigs_10k.fasta \
        --coverage_file {output.cov_table} \
        -b {config[outdir]}/{wildcards.sample}/binning/concoct.out/results 2> {log}

        merge_cutup_clustering.py \
        {config[outdir]}/{wildcards.sample}/binning/concoct.out/results/clustering_gt1000.csv \
        > {output.merged_csv}

        mkdir -p {config[outdir]}/{wildcards.sample}/binning/concoct.out/results/fasta_bins
        extract_fasta_bins.py {input.contigs_filt} {output.merged_csv} \
        --output_path {config[outdir]}/{wildcards.sample}/binning/concoct.out/results/fasta_bins
        """

rule maxbin:
    input:
        hr1 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_1_hr.fastq.gz"),
        hr2 = os.path.join(config["outdir"], "{sample}", "preprocessing", "{sample}_2_hr.fastq.gz"),
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta")
    threads: 24
    conda: "../envs/maxbin_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "binning", "maxbin.out"))
    log:
        os.path.join(config["outdir"], "logs", "maxbin", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "maxbin", "{sample}_bmrk.txt")
    shell:
        """
        
        run_MaxBin.pl -thread {threads} -contig {input.contigs_filt} \
        -reads {input.hr1} -reads2 {input.hr2} \
        -out {config[outdir]}/{wildcards.sample}/binning/maxbin.output 2> {log} || true

        mkdir -p {output}
        mv {config[outdir]}/{wildcards.sample}/binning/maxbin.output* {output}
        """

rule metabat:
    input:
        contigs_filt = os.path.join(config["outdir"], "{sample}", "assembly", "contigs_filt_1000bp.fasta"),
        sorted_bam = os.path.join(config["outdir"], "{sample}", "binning", "map_reads", "sorted_bam", "{sample}_sorted.bam")
    threads: 24
    conda: "../envs/metabat_env.yaml"
    output:
        depth = os.path.join(config["outdir"], "{sample}", "binning", "depth.txt"),
        outfile = directory(os.path.join(config["outdir"], "{sample}", "binning", "metabat.out"))
    log:
        os.path.join(config["outdir"], "logs", "metabat", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "metabat", "{sample}_bmrk.txt")
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.sorted_bam}
        mkdir -p {output.outfile}
        
        metabat2 -i {input.contigs_filt} -a {output.depth} -o {output.outfile}/bin -v 2> {log}
        """

rule tsv_files_for_dastool:
    input:
        concoct = os.path.join(config["outdir"], "{sample}", "binning", "concoct.out", "results", "clustering_merged.csv"),
        maxbin = os.path.join(config["outdir"], "{sample}", "binning", "maxbin.out"),
        metabat = os.path.join(config["outdir"], "{sample}", "binning", "metabat.out")
    conda: "../envs/dastool_env.yaml"
    output:
        concoct_tsv = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "concoct.contigs2bin.tsv"),
        maxbin_tsv = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "maxbin.contigs2bin.tsv"),
        metabat_tsv = os.path.join(config["outdir"], "{sample}", "binning", "dastool", "metabat.contigs2bin.tsv")
    shell:
        """
        mkdir -p {config[outdir]}/{wildcards.sample}/binning/dastool

        # Convert concoct csv to tsv
        perl -pe 's/,/\tCONCOCT.bin./g' {input.concoct} > {output.concoct_tsv}
        sed -i '/^contig/d' {output.concoct_tsv}

        # Convert maxbin bins to tsv
        Fasta_to_Contig2Bin.sh -i {input.maxbin} -e fasta > {output.maxbin_tsv}

        # Convert metabat bins to tsv
        Fasta_to_Contig2Bin.sh -i {input.metabat} -e fa > {output.metabat_tsv}
        """
        
