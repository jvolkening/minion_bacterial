Companion materials for the manuscript "Taylor et al: *Rapid, multiplexed, whole genome and plasmid sequencing of foodborne pathogens using long-read nanopore technology*"
=======================================================

This repository contains the primary Nextflow workflow (`minion_bacterial.nf`)
for bacterial nanopore assembly described in the manuscript and accessory
scripts for reproducing the results therein. Workflows for the phylogenetic
analysis can be run separately from the `phylogenetics` directory. The time
series analysis can be reproduced using the following steps:


1. Initialize the software environments used in the build script. A separate
environment is needed for Poretools as it's conda package is incompatible with
Porechop.

```
conda env create -n minion_init -f conda_init.yml
conda create --name poretools poretools
```

2. Download the raw data, extract and trim the FASTQ reads, and create the
time series subsampling directories for use in benchmarking different MinION
run times:

```
./init.sh
```

3. Run the run-time benchmarking on the *Salmonella* and *E. coli* datasets

```
bin/bm_Sent.pl --in ts/Sent --out <path_to_output_dir> --meta meta --threads <num_threads>
bin/bm_Ecoli.pl --in ts/Ecoli --out <path_to_output_dir> --meta meta --threads <num_threads>
```

## NOTES

### Computing resources

The analyses for the manuscript were run on a server with 64 cores and 512 GB
memory. Certain steps (`nanopolish` in particular) are quite slow even with 64
cores available (see manuscript for run times and other benchmarking details).
Be aware of this if attempting to replicate the analysis on low-end hardware.

### Read trimming

Please note that the version of Porechop (0.2.2) used in the manuscript to
demultiplex and trim the reads is older than the earliest version available in
Bioconda (0.2.3). It is possible that this may result in slight differences in
the output metrics when reproducing the manuscript analysis.

### Unicycler versions

If you generate assemblies with significantly poorer QC metrics than those
described in the manuscript, please check the Unicycler logs in your Nextflow
results directory. All of the conda software packages in the workflow are
versioned identically to those used in the manuscript. However, we experienced
issues with Unicycler when the associated `racon` dependency was silently
upgraded from v1.3.1 to v1.3.3 -- specifically, the same version of
Unicycler would now silently fail to perform Racon polishing and would output
the unpolished assembly. You can see this in the Unicycler log, as no rounds
of Racon polishing will be reported in that section. Pinning the conda Racon
version to 1.3.1 (the version used for the manuscript) in the Nextflow script
seems to have fixed this issue, but it is something to keep an eye on if you
attempt to use the workflow and experience poor results.

### Nanopolish and FAST5 file paths

Some tools produce FASTQ files with the path to the FAST5 files from which
each read was derived embedded in the read descriptor. The version of
Nanopolish used in this manuscript appears to use that information in
preference to the summary table provided on the command line. If those FAST5
files cannot be found at the same path (absolute or relative) found in the
FASTQ headers, Nanopolish will give a series of warnings but will not die.
Rather, it will output a "polished" assembly identical to the input. Since
Nanopolish is being run within a Nextflow workflow, those warnings are easily
missed.

In the `init.sh` script used here, the absolute path to the FAST5 directory is
provided to `poretools` to ensure that the same absolute path is encoded in
the FASTQ headers. If you are using this workflow on your own data, please
double-check that these paths are correctly encoded for Nanopolish to find
them. We have found this to be tricky to troubleshoot, but the most obvious
sign that something is wrong is that the Nanopolish input and output
assemblies will be identical.

### SRA FAST5 data

The SRA accessions for this project are SRR9603470 (Salmonella) and SRR9603471 (E. coli). The FAST5 data can be downloaded via the SRA web portal under the 'Download' tab. However, be aware that, despite the file naming SRA used, these files are tarballs only and not gzipped, so set your `tar` extraction flags appropriately. You can find example commands for fetching the data in the `init.sh` script in this repo.
