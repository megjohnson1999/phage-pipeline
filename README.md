## How to run:


snakemake \
--profile ../profile/slurm/ \
--config [options]


## options:

reads: specify path to directory where paired-end fastq reads are

outdir: specify path to directory where all outputs will be created

fastp_min_sequence_length: length threshold (in bp) for fastp step (default is 120)

human_ref: "/ref/sahlab/data/GRCh38.fna.gz"

genomad_database: "/ref/sahlab/data/viral_analysis_DBs/genomad_DBs/genomad_db"

bakta_database: "/ref/sahlab/data/bakta_db"

cat_database: "/ref/sahlab/data/CAT_prepare_20210107"


## Example:

snakemake --profile ../profile/slurm/ --config reads=/scratch/sahlab/Megan/test_reads outdir=/scratch/sahlab/Megan/pipeline_test_out
