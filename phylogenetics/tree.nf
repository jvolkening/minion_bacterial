#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.minion         = null
params.illumina       = null
params.reference      = null
params.threads        = 1
params.outdir         = 'results'

reference = file(params.reference)

minion_vcf = Channel.fromPath( params.minion )
illumina_vcf = Channel.fromPath( params.illumina )


/* ***************************************************************************
   STEP 1a :: Illumina filtering
*************************************************************************** */
process filter_illumina {

    tag "$vcf"

    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    file vcf from illumina_vcf

    output:
    file "*.filt.vcf" into illumina_filt

    script:
    """
    filter_vcf.pl \
        --in $vcf \
        --out ${vcf.baseName}.filt.vcf \
        --drop_indels
    """

}

/* ***************************************************************************
   STEP 1b :: MinION tagging
*************************************************************************** */
process tag_minion {

    tag "$vcf"

    publishDir "${params.outdir}/tagged", mode: 'copy'

    conda 'bioconda::perl-biox-seq'

    input:
    each file(vcf) from minion_vcf
    file ref from reference

    output:
    file "*.filt.vcf" into minion_filt

    script:
    """
    tag_snps.pl \
        --ref $ref \
        --in $vcf \
        --out ${vcf.baseName}.filt.vcf  \
        --tag_agtc
    """

}

/* ***************************************************************************
   STEP 2 :: Generate SNP matrix
*************************************************************************** */
process make_matrix {

    tag "$illumina_vcfs AND $minion_vcfs"

    conda 'bioconda::perl-biox-seq'

    publishDir "${params.outdir}/trees", mode: 'copy'

    input:
    file minion_vcfs from minion_filt.collect()
    file illumina_vcfs from illumina_filt.collect()
    file ref from reference

    output:
    file "matrix.fa" into snp_matrix
    file "matrix.log"

    script:
    """
    gen_matrix.pl \
        --ref $ref \
        $minion_vcfs \
        $illumina_vcfs \
    > matrix.fa \
    2> matrix.log
    """

}

/* ***************************************************************************
   STEP 2 :: Generate ML tree with SH support
*************************************************************************** */
process phyml_sh {

    tag "$matrix"

    conda 'bioconda::phyml=3.3.20190321 bioconda::perl-biox-seq'

    publishDir "${params.outdir}/tree", mode: 'copy'

    input:
    file matrix from snp_matrix

    output:
    file "sh.newick"

    script:
    """
    fasta2phylip < $matrix > aln.phylip
    phyml -i aln.phylip -b -4
    mv aln.phylip_phyml_tree.txt sh.newick
    """

}

/* ***************************************************************************
   STEP 2 :: Generate ML tree with bootstrap support
*************************************************************************** */
process phyml_bs {

    tag "$matrix"

    conda 'bioconda::phyml=3.3.20190321 bioconda::perl-biox-seq'

    publishDir "${params.outdir}/tree", mode: 'copy'

    input:
    file matrix from snp_matrix

    output:
    file "bs500.newick"

    script:
    """
    fasta2phylip < $matrix > aln.phylip
    phyml -i aln.phylip -b 500
    #mpiexec phyml -i aln.phylip -b 500
    mv aln.phylip_phyml_tree.txt bs500.newick
    """

}
