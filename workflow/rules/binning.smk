rule binning_prep:
    input:
        hr1 = "reads/host_removed/{sample}_1_hr.fastq.gz",
        hr2 = "reads/host_removed/{sample}_2_hr.fastq.gz",
        contigs = "out/{sample}/assembly/contigs.fasta"
    threads: 8
    conda: "../envs/bowtie2.yaml"
    output:
        contigs_filt = "out/{sample}/assembly/contigs_filt_1000bp.fasta",
        sam = "out/{sample}/binning/map_reads/{sample}.sam",
        bam = "out/{sample}/binning/map_reads/{sample}.bam",
        sorted_bam = "out/{sample}/binning/map_reads/sorted_bam/{sample}_sorted.bam"
    log:
        "logs/binning_prep/{sample}.log"
    shell:
        """
        cat {input.contigs} | seqkit seq -m 1000 > {output.contigs_filt}

        bowtie2-build {output.contigs_filt} out/{wildcards.sample}/binning/map_reads/btdb &>> {log}

        bowtie2 --no-unal -p {threads} -x out/{wildcards.sample}/binning/map_reads/btdb -1 {input.hr1} -2 {input.hr2} -S {output.sam} &>> {log}

        samtools view -@ {threads} -Sb -o {output.bam} {output.sam} &>> {log}

        samtools sort -O bam -o {output.sorted_bam} {output.bam} &>> {log}

        samtools index {output.sorted_bam} &>> {log}
        """

rule concoct:
    input:
        contigs_filt = "out/{sample}/assembly/contigs_filt_1000bp.fasta",
        sorted_bam = "out/{sample}/binning/map_reads/sorted_bam/{sample}_sorted.bam"
    threads: config.get("num_threads", 8)
    conda: "../envs/concoct_env.yaml"
    output:
        bed = "out/{sample}/binning/concoct.out/contigs_10k.bed",
        fa = "out/{sample}/binning/concoct.out/contigs_10k.fasta",
        cov_table = "out/{sample}/binning/concoct.out/coverage_table.tsv",
        merged_csv = "out/{sample}/binning/concoct.out/results/clustering_merged.csv"
    log:
        "logs/concoct/{sample}.log"
    shell:
        """
        mkdir -p out/{wildcards.sample}/binning/concoct.out
        mkdir -p out/{wildcards.sample}/binning/concoct.out/results

        cut_up_fasta.py {input.contigs_filt} -c 10000 -o 0 --merge_last -b {output.bed} > {output.fa} 2>> {log}

        concoct_coverage_table.py out/{wildcards.sample}/binning/concoct.out/contigs_10k.bed {input.sorted_bam} > {output.cov_table} 2>> {log}

        concoct --composition_file out/{wildcards.sample}/binning/concoct.out/contigs_10k.fasta --coverage_file {output.cov_table} \
        -b out/{wildcards.sample}/binning/concoct.out/results 2>> {log}

        merge_cutup_clustering.py out/{wildcards.sample}/binning/concoct.out/results/clustering_gt1000.csv > {output.merged_csv} 2>> {log}

        mkdir -p out/{wildcards.sample}/binning/concoct.out/results/fasta_bins
        extract_fasta_bins.py {input.contigs_filt} {output.merged_csv} --output_path out/{wildcards.sample}/binning/concoct.out/results/fasta_bins 2>> {log}
        """

rule maxbin:
    input:
        hr1 = "reads/host_removed/{sample}_1_hr.fastq.gz",
        hr2 = "reads/host_removed/{sample}_2_hr.fastq.gz",
        contigs_filt = "out/{sample}/assembly/contigs_filt_1000bp.fasta"
    threads: config.get("num_threads", 8)
    conda: "../envs/maxbin_env.yaml"
    output:
        directory("out/{sample}/binning/maxbin.out")
    log:
        "logs/maxbin/{sample}.log"
    shell:
        """
        run_MaxBin.pl -thread {threads} -contig {input.contigs_filt} \
        -reads {input.hr1} -reads2 {input.hr2} \
        -out out/{wildcards.sample}/binning/maxbin.output 2>> {log} || true

        mkdir -p {output}
        mv out/{wildcards.sample}/binning/maxbin.output* {output} 2>> {log}
        """

rule metabat:
    input:
        contigs_filt = "out/{sample}/assembly/contigs_filt_1000bp.fasta",
        sorted_bam = "out/{sample}/binning/map_reads/sorted_bam/{sample}_sorted.bam"
    threads: config.get("num_threads", 8)
    conda: "../envs/metabat_env.yaml"
    output:
        depth = "out/{sample}/binning/depth.txt",
        outfile = directory("out/{sample}/binning/metabat.out")
    log:
        "logs/metabat/{sample}.log"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.sorted_bam}
        mkdir -p {output.outfile}
        #cd {output.outfile}
        #runMetaBat.sh ../../../{input.contigs_filt} ../../../{input.sorted_bam} 2> {log}
        metabat2 -i {input.contigs_filt} -a {output.depth} -o {output.outfile}/bin -v 2>> {log}
        """

rule tsv_files_for_dastool:
    input:
        concoct = "out/{sample}/binning/concoct.out/results/clustering_merged.csv",
        maxbin = "out/{sample}/binning/maxbin.out",
        metabat = "out/{sample}/binning/metabat.out"
    conda: "../envs/dastool_test2.yaml"
    output:
        concoct_tsv = "out/{sample}/binning/dastool/concoct.contigs2bin.tsv",
        maxbin_tsv = "out/{sample}/binning/dastool/maxbin.contigs2bin.tsv",
        metabat_tsv = "out/{sample}/binning/dastool/metabat.contigs2bin.tsv"
    shell:
        """
        mkdir -p "out/{wildcards.sample}/binning/dastool"

        # Convert concoct csv to tsv
        perl -pe 's/,/\tCONCOCT.bin./g' {input.concoct} > {output.concoct_tsv}
        sed -i '/^contig/d' {output.concoct_tsv}

        # Convert maxbin bins to tsv
        Fasta_to_Contig2Bin.sh -i {input.maxbin} -e fasta > {output.maxbin_tsv}

        # Convert metabat bins to tsv
        Fasta_to_Contig2Bin.sh -i {input.metabat} -e fa > {output.metabat_tsv}
        """

