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

if (params.illumina != null) {
    Channel
        .fromFilePairs( params.illumina, size: 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching" }
        .set{ to_trim }
}

if (params.minion != null) {
    to_sim = Channel.fromPath( params.minion )
}


/* ***************************************************************************
   STEP 1a :: Read trimming (Illumina)
*************************************************************************** */
if (params.illumina != null) {
    process trim {

        tag "$reads"

        conda 'bioconda::trim-galore=0.5.0'

        publishDir "${params.outdir}/trim", mode: 'copy'

        input:
        set val(id),file(reads) from to_trim

        output:
        set val(id),file("*{1,2}.fq.gz") into to_align

        script:
        """
        trim_galore \
            --gzip \
            --paired \
            --trim-n \
            --quality 15 \
            --stringency 5 \
            --length 30 \
            --max_n 4 \
            $reads
        """

    }
}

/* ***************************************************************************
   STEP 1b :: Read simulation (MinION)
*************************************************************************** */
if (params.minion != null) {
    process simulate {

        tag "$assembly"

        conda 'bioconda::art=2016.06.05'

        //publishDir "${params.outdir}/trim", mode: 'copy'

        input:
        file assembly from to_sim

        output:
        set val(assembly.baseName),file("*{1,2}.fq.gz") into to_align

        script:
        """
        art_illumina \
            -i $assembly \
            --fcov 50 \
            -ss MSv3 \
            -l 150 \
            -p \
            -m 300 \
            -s 50 \
            -na \
            -o ${assembly.baseName}_
        gzip *.fq
        """

    }
}

/* ***************************************************************************
   STEP 2 :: Read mapping
*************************************************************************** */
process align {

    tag "$reads on $ref"
    cpus "${params.threads}"
    memory '12 G'
 
    conda 'bioconda::bwa=0.7.17 ncurses bioconda::samtools=1.9'

    publishDir "${params.outdir}/aligned", mode: 'copy'

    input:
    set val(id),file(reads) from to_align
    file ref from reference

    output:
    set val(id),file("${id}.bam") into to_call

    script:
    """
    bwa index $ref
    bwa mem \
        -t ${params.threads} \
        $ref \
        $reads \
    | samtools view -Sbh -F 4 \
    | samtools sort \
    > ${id}.bam
    """

}

/* ***************************************************************************
   STEP 3 :: Variant calling
*************************************************************************** */
process call {

    tag "$bam on $ref"
    cpus "${params.threads}"
    memory '4 G'

    conda 'ncurses bioconda::samtools=1.9 bioconda::lofreq=2.1.3.1'

    publishDir "${params.outdir}/called", mode: 'copy'

    input:
    set val(id),file(bam) from to_call
    file ref from reference

    output:
    file "${id}.vcf"

    script:
    """
    samtools index $bam
    lofreq indelqual \
        --dindel \
        -f $ref \
        -o dindel.bam \
        $bam
    samtools index dindel.bam
    lofreq call-parallel \
        --call-indels \
        --pp-threads ${params.threads} \
        -f $ref \
        -o ${id}.vcf \
        dindel.bam
    """

}
