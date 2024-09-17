
![updated_workflow_091024 drawio](https://github.com/user-attachments/assets/4f59c6d9-a453-4985-a11e-f8eed6714539)


# How to run:
Set up an environment with Snakemake version 8+, [mamba](https://anaconda.org/conda-forge/mamba), and [snakemake-executor-plugin-slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)

```
cd workflow

snakemake --profile ../profile/slurm/ --config [options]
```

# Options:

 - reads: specify path to directory where paired-end fastq reads are (with suffixes _1.fastq.gz and _2.fastq.gz)

 - outdir: specify path to directory where all outputs will be created

 - fastq_names_1: default is {sample}_1.fastq.gz

 - fastq_names_2: default is {sample}_2.fastq.gz

 - fastp_min_sequence_length: length threshold (in bp) for fastp step (default is 120)

 - human_ref: "/ref/sahlab/data/GRCh38.fna.gz"

 - genomad_database: "/ref/sahlab/data/viral_analysis_DBs/genomad_DBs/genomad_db"

 - bakta_database: "/ref/sahlab/data/bakta_db"

 - cat_database: "/ref/sahlab/data/CAT_prepare_20210107"

 - checkv_database: "/ref/sahlab/data/viral_analysis_DBs/checkV_DB/checkv-db-v1.4"


### Example command:

```
snakemake --profile ../profile/slurm/ --config reads=/scratch/sahlab/Megan/test_reads outdir=/scratch/sahlab/Megan/pipeline_test_out
```

# Outputs:

The output directory should contain separate directories for each sample. Each sample's directory should have 5 subdirectories:

 - assembly

 - binning

 - coverm

 - taxonomy

 - phage_analysis

### In the phage_analysis directory:

 - final_prophage_table.tsv: has the prophages (contigs and start/stop coordinates)

 - final_prophage_table_with_taxonomy.tsv: has the prophages + taxonomy info (contigs, start/stop coordinates, and taxonomy assigned to that contig)

 - final_prophage.fasta: has the sequences of all the final prophage regions

