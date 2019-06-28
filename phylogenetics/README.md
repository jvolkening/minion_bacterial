This subdirectory contains code to reproduce the *Salmonella* phylogenetic
analyses presented in the manuscript. To replicate the results, please
follow these steps:

1. Clone the reference dataset repository from the benchmark paper
(<https://dx.doi.org/10.7717%2Fpeerj.3893>)

```
git clone https://github.com/WGS-standards-and-analysis/datasets
```

2. Within that repo, use their script to download the reference Illumina
datasets. We need to use a patched data table with updated accessions and
checksums.

```
cd datasets
scripts/GenFSGopher.pl \
    --outdir test_data \
    ../misc/Salmonella_enterica_1203NYJAP-1.patched.tsv
cd ..
```

3. Copy the MinION assembly FASTA files to analyze into the `minion_assemblies`
directory.

4. Run variant calling for MinION assemblies and Illumina reference data,
setting the number of available threads as appropriate.

```
nextflow run snps.nf \
    --threads <num_threads> \
    --reference datasets/test_data/CFSAN000212.fasta \
    --minion 'minion_assemblies/*.fa' \
    --outdir minion_res
nextflow run snps.nf \
    --threads <num_threads> \
    --reference datasets/test_data/CFSAN000212.fasta \
    --illumina 'datasets/test_data/*_{1,2}.fastq.gz' \
    --outdir illumina_res
```

5. Generate tree.

```
nextflow run tree.nf \
    --threads <num_threads> \
    --reference datasets/test_data/CFSAN000212.fasta \
    --minion 'minion_res/call/*.vcf' \
    --illumina 'illumina_res/call/*.vcf' \
    --outdir tree_res
```

Trees will be in the output directory (here, `tree_res/trees`).
