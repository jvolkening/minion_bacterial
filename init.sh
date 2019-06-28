#!/usr/bin/env bash

# Fetch BUSCO dataset
wget -P meta http://busco.ezlab.org/v2/datasets/enterobacteriales_odb9.tar.gz
tar -C meta -xvf meta/enterobacteriales_odb9.tar.gz

source activate minion_init

# Fetch SRA Illumina data
fastq-dump --gzip --split-files --skip-technical --clip --dumpbase --outdir raw/illumina/Ecoli SRR6373397

# Fetch SRA MinION data
## fetch BC01 to raw/minion/Sent/fast5
## fetch BC02 to raw/minion/Ecoli/fast5

THREADS=1
EC_UNTRIMMED=raw/minion/Ecoli/untrimmed.fq.gz
EC_TRIMMED=raw/minion/Ecoli/trimmed.fq.gz
SE_UNTRIMMED=raw/minion/Sent/untrimmed.fq.gz
SE_TRIMMED=raw/minion/Sent/trimmed.fq.gz

# Extract FASTQ (the readlink is added so that absolute paths are used in the
# FASTQ header. Nanopolish seems to use this information during indexing.

source deactivate
source activate poretools

poretools fastq $(readlink -f raw/minion/Sent/fast5)  | gzip > $SE_UNTRIMMED
poretools fastq $(readlink -f raw/minion/Ecoli/fast5) | gzip > $EC_UNTRIMMED

source deactivate
source activate minion_init

# Trim adapters (manuscript used porechop v0.2.2, which is not available in Bioconda
# Note that a few reads (~1%) will be slightly different at the ends between
# those output below and those used in the manuscript. I believe this is
# because originally the reads were both demultiplexed and trimmed by
# Porechop, and thus Porechop searched for both barcodes in all reads. Because
# they needed to be demultiplexed as FAST5 before submission to SRA, this
# difference was unavoidable.
porechop --input $SE_UNTRIMMED \
    --threads $THREADS \
    --adapter_threshold 95.0 \
    --middle_threshold 100 \
    --check_reads 4000 \
| gzip > $SE_TRIMMED
porechop --input $EC_UNTRIMMED \
    --threads $THREADS \
    --adapter_threshold 95.0 \
    --middle_threshold 100 \
    --check_reads 4000 \
| gzip > $EC_TRIMMED


# Create subsampled datasets based on elapsed time series
bin/time_series.pl \
    --read_map meta/Sent.map.tsv.gz \
    --in_fast5 raw/minion/Sent/fast5 \
    --in_trim $SE_TRIMMED \
    --in_untrim $SE_UNTRIMMED \
    --durations 15,30,60,120,240,480,960,1500 \
    --out_dir ts/Sent
bin/time_series.pl \
    --read_map meta/Ecoli.map.tsv.gz \
    --in_fast5 raw/minion/Ecoli/fast5 \
    --in_trim $EC_TRIMMED \
    --in_untrim $EC_UNTRIMMED \
    --durations 15,30,60,120,240,480,960,1500 \
    --out_dir ts/Ecoli
