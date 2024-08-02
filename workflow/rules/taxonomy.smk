rule cat:
    input:
        os.path.join(config["outdir"], "{sample}", "binning", "final_filtered_contigs.fasta")
    params:
        db = os.path.join(config["cat_database"], "2021-01-07_CAT_database"),
        tax = os.path.join(config["cat_database"], "2021-01-07_taxonomy")
    threads: 24
    conda: "../envs/cat_env.yaml"
    output:
        directory(os.path.join(config["outdir"], "{sample}", "taxonomy", "CAT"))
    log:
        os.path.join(config["outdir"], "logs", "CAT", "{sample}.log")
    benchmark:
        os.path.join(config["outdir"], "benchmarks", "CAT", "{sample}_bmrk.txt")
    shell:
        """
        set -ue
        mkdir -p {output}

        # Main CAT pipeline
        CAT contigs -c {input} -o {output}/out.CAT \
        --force --sensitive \
        -d {params.db} -t {params.tax} -n {threads}

        # Add taxonomic names
        CAT add_names -i {output}/out.CAT.contig2classification.txt \
        -o {output}/out.CAT.taxonomy \
        -t {params.tax} --only_official --exclude_scores --force

        # Summarize results
        CAT summarise -c {input} -i {output}/out.CAT.taxonomy \
        -o {output}/out.CAT.summary

        # Create contig taxonomy table
        cut -f1,6-12 {output}/out.CAT.taxonomy > {output}/contig.taxonomy
        cut -f1-5 {output}/out.CAT.taxonomy > {output}/CAT.scores
        """
