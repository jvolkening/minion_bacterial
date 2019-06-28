#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.trimmed        = null
params.untrimmed      = null
params.indexed        = null
params.outdir         = null
params.reference      = null
params.fast5dir       = null
params.buscolineage   = null
params.illumina       = null
params.id_map         = null
params.genomelength   = 5000000
params.depth          = 50
params.threads        = 1
params.qweight        = 5
params.run_pilon      = false
params.starts         = null
params.minimus2       = "${baseDir}/bin/minimus2_fixed"


buscolineage = file(params.buscolineage)
minimus2     = file(params.minimus2)
trimmed      = file(params.trimmed)
if (params.untrimmed) {
    untrimmed = file(params.untrimmed)
}
if (params.starts) {
    starts = file(params.starts)
}
if (params.id_map) {
    id_map    = file(params.id_map)
}
if (params.reference) {
    reference = file(params.reference)
}

if (params.illumina) {
    Channel.fromPath(params.illumina).into { illumina_lofreq; illumina_pilon }
}


/* ***************************************************************************
   STEP 1 :: Quality/length filtering
*************************************************************************** */

process filtlong {

    tag "$reads"

    conda 'bioconda::filtlong=0.2.0'

    publishDir "${params.outdir}/filtlong", mode: 'copy'

    input:
    file reads from trimmed

    output:
    file "${reads.baseName}.filt.fq" into filtered

    script:
    def depth = params.genomelength * params.depth
    """
    filtlong -t $depth --mean_q_weight ${params.qweight} $reads \
      > ${reads.baseName}.filt.fq
    """

}

/* ***************************************************************************
   STEP 2 :: Assembly
*************************************************************************** */

process unicycler {

    tag "$reads"

    conda 'bioconda::unicycler=0.4.6 bioconda::racon=1.3.1'

    publishDir "${params.outdir}/unicycler", mode: 'copy'

    input:
    file reads from filtered

    output:
    file "assembly.gfa" into assembly
    file "assembly.fasta"
    file "unicycler.log"

    script:
    """
    unicycler -l $reads -o . --keep 0 -t ${params.threads} --no_rotate
    """

}

/* ***************************************************************************
   STEP 3 :: Circularize
*************************************************************************** */

process circularize {

    tag "$assembly"

    conda 'bioconda::amos bioconda::perl-biox-seq'

    publishDir "${params.outdir}/circularized", mode: 'copy'

    input:
    file assembly from assembly
    file minimus2 from minimus2

    output:
    file "circularized.fasta" into circularized
    file "circularized.log"

    script:
    """
    circularize.pl $minimus2 < $assembly 2> circularized.log > circularized.fasta
    """

}

if (!params.starts) {
    final_assembly = circularized
} else {

    /* ***************************************************************************
    STEP 4 :: Orient
    *************************************************************************** */

    process orient {

        tag "$assembly"

        conda 'bioconda::bwa=0.7.17 bioconda::perl-biox-seq'

        publishDir "${params.outdir}/orient", mode: 'copy'

        input:
        file assembly from circularized

        output:
        file "oriented.fasta" into final_assembly
        file "oriented.log"

        script:
        """
        orient.pl --assembly $assembly --starts $starts \
          --threads ${params.threads} \
          2> oriented.log \
          > oriented.fasta
        """

    }
}


if (params.fast5dir) {

    /* ***************************************************************************
    STEP 4A :: Nanopolish
    *************************************************************************** */
    
    process preprocess {

        tag "$ref"

        conda 'ncurses bioconda::nanopolish=0.10.2 bioconda::minimap2=2.12 bioconda::samtools=1.9'

        input:
        file ref from final_assembly
        file reads from untrimmed

        output:
        file ref into np_ref, vcf_ref
        file 'untrimmed.bam' into np_bam
        file 'untrimmed.bam.bai' into np_bam_index
        file reads into np_reads
        file "${reads}.*" into np_index
        file 'ranges.txt' into np_ranges

        script:
        def indexed = params.indexed ?: params.untrimmed
        def n_thr = params.threads > 4 ? 4 : params.threads
        """
        minimap2 -d ${ref}.mmi $ref
        minimap2 -ax map-ont -t ${params.threads} ${ref}.mmi  $reads \
        | samtools sort -@${params.threads} -o untrimmed.bam
        samtools index untrimmed.bam
        zcat $id_map | cut -f1,2 > summary
        nanopolish index -d ${params.fast5dir} $reads \
          --sequencing-summary summary
        nanop=\$(type -P nanopolish)
        python \$(dirname \$nanop)/nanopolish_makerange.py $ref > ranges.txt
        """

    }

    np_ranges = np_ranges.splitText().map{it -> it.trim()}

    process nanopolish {

        tag "$ref $win"
        maxForks = params.threads > 4 ? (int) (params.threads/4) : 1
        // nanopolish sometimes fails inexplicably, but works when re-tried
        errorStrategy 'retry'
        maxErrors 2

        conda 'bioconda::nanopolish=0.10.2'

        publishDir "${params.outdir}/nanopolish/vcf", mode: 'copy'

        input:
        file reads from np_reads
        file index from np_index.collect()
        val win from np_ranges
        file ref from np_ref
        file bam from np_bam
        file bam_index from np_bam_index

        output:
        file {"${win}.vcf".replaceAll(/:/, ".")} into np_vcf

        script:
        def n_thr = params.threads > 4 ? 4 : params.threads
        def vcf = "${win}.vcf".replaceAll(/:/, ".")
        """
        nanopolish variants -o ${vcf} -w $win -r $reads \
          -b $bam -g $ref --threads $n_thr \
          --consensus \
          --methylation-aware dcm,dam \
          --fix-homopolymers
        """

    }

    process vcf2fasta {

        tag "$ref"

        conda 'bioconda::nanopolish=0.10.2'

        publishDir "${params.outdir}/nanopolish", mode: 'copy'

        input:
        file vcfs from np_vcf.collect()
        file ref from vcf_ref

        output:
        file "${ref.baseName}.polished.fa" into np_2
        file "${ref.baseName}.polished.fa" into polished_assembly_1

        script:
        """
        nanopolish vcf2fasta -g $ref $vcfs > ${ref.baseName}.polished.fa
        """

    }

    final_assembly = final_assembly.mix(polished_assembly_1)

    process preprocess_2 {

        tag "$ref"

        conda 'ncurses bioconda::nanopolish=0.10.2 bioconda::minimap2=2.12 bioconda::samtools=1.9'

        input:
        file ref from np_2
        file reads from untrimmed

        output:
        file ref into np_ref_2, vcf_ref_2
        file 'untrimmed.bam' into np_bam_2
        file 'untrimmed.bam.bai' into np_bam_index_2
        file reads into np_reads_2
        file "${reads}.*" into np_index_2
        file 'ranges.txt' into np_ranges_2

        script:
        def indexed = params.indexed ?: params.untrimmed
        def n_thr = params.threads > 4 ? 4 : params.threads
        """
        minimap2 -d ${ref}.mmi $ref
        minimap2 -ax map-ont -t ${params.threads} ${ref}.mmi  $reads \
        | samtools sort -@${params.threads} -o untrimmed.bam
        samtools index untrimmed.bam
        zcat $id_map | cut -f1,2 > summary
        nanopolish index -d ${params.fast5dir} $reads \
          --sequencing-summary summary
        nanop=\$(type -P nanopolish)
        python \$(dirname \$nanop)/nanopolish_makerange.py $ref > ranges.txt
        """

    }

    np_ranges_2 = np_ranges_2.splitText().map{it -> it.trim()}

    process nanopolish_2 {

        tag "$ref $win"
        maxForks = params.threads > 4 ? (int) (params.threads/4) : 1
        // nanopolish sometimes fails inexplicably, but works when re-tried
        errorStrategy 'retry'
        maxErrors 2

        conda 'nanopolish=0.10.2'

        publishDir "${params.outdir}/nanopolish/vcf", mode: 'copy'

        input:
        file reads from np_reads_2
        file index from np_index_2.collect()
        val win from np_ranges_2
        file ref from np_ref_2
        file bam from np_bam_2
        file bam_index from np_bam_index_2

        output:
        file {"${win}.vcf".replaceAll(/:/, ".")} into np_vcf_2

        script:
        def n_thr = params.threads > 4 ? 4 : params.threads
        def vcf = "${win}.vcf".replaceAll(/:/, ".")
        """
        nanopolish variants -o ${vcf} -w $win -r $reads \
          -b $bam -g $ref --threads $n_thr \
          --consensus \
          --methylation-aware dcm,dam \
          --fix-homopolymers
        """

    }

    process vcf2fasta_2 {

        tag "$ref"

        conda 'nanopolish=0.10.2'

        publishDir "${params.outdir}/nanopolish", mode: 'copy'

        input:
        file vcfs from np_vcf_2.collect()
        file ref from vcf_ref_2

        output:
        file "${ref.baseName}.polished.fa" into polished_assembly_2

        script:
        """
        nanopolish vcf2fasta -g $ref $vcfs > ${ref.baseName}.polished.fa
        """

    }

    final_assembly = final_assembly.mix(polished_assembly_2)

}
else if (params.run_pilon) {

    /* ***************************************************************************
    STEP 4B :: pilon
    *************************************************************************** */

    process pilon {

        tag "$ref"

        conda 'bioconda::bwa=0.7.17 ncurses bioconda::samtools=1.9 bioconda::pilon=1.2.2 bioconda::perl-biox-seq'

        publishDir "${params.outdir}/pilon", mode: 'copy'

        input:
        file ref from final_assembly
        file reads from illumina_pilon.collect()

        output:
        file "${ref.baseName}.pilon.fasta" into polished_assembly

        script:
        """
        bwa index $ref
        bwa mem -t ${params.threads} $ref $reads \
            | samtools sort -@${params.threads} -o sorted.bam -
        samtools index sorted.bam
        pilon -Xmx100G --genome $ref --bam sorted.bam --fix snps,indels \
          --output pilon
        copy_headers.pl $ref pilon.fasta > ${ref.baseName}.polished.fasta
        """
    }

    final_assembly = final_assembly.mix(polished_assembly)

}


final_assembly.into {
    final_1
    final_2
}

/* ***************************************************************************
STEP 5 :: BUSCO
*************************************************************************** */

process busco {

    tag "$assembly"

    conda 'bioconda::busco=3.0.2 bioconda::blast=2.7.1'

    publishDir "${params.outdir}/busco", mode: 'copy'
    afterScript "for i in run_${assembly.baseName}/*; do mv \$i .; done; rmdir run_${assembly.baseName}"

    input:
    file assembly from final_2

    output:
    file "full_table_${assembly.baseName}.tsv"

    script:
    """
    run_busco -i $assembly -o ${assembly.baseName} -m genome -l $buscolineage \
      -c ${params.threads} --blast_single_core
    """
}


/* ***************************************************************************
STEP 6 :: dnadiff
*************************************************************************** */

if (params.reference) {
    process dnadiff {

        tag "$assembly"

        conda 'bioconda::mummer=3.23'

        publishDir "${params.outdir}/dnadiff", mode: 'copy'

        input:
        file assembly from final_1
        file ref from reference

        output:
        file "${assembly.baseName}.1coords"
        file "${assembly.baseName}.report"

        script:
        """
        dnadiff -p ${assembly.baseName} $ref $assembly 
        """
    }
}

/* ***************************************************************************
STEP 7 :: lofreq
*************************************************************************** */

else if (params.illumina) {
    process lofreq {

        tag "$assembly with $reads"

        conda 'bioconda::bwa=0.7.17 ncurses bioconda::samtools=1.9 bioconda::bedtools=2.27.1 bioconda::lofreq=2.1.3.1'

        publishDir "${params.outdir}/lofreq", mode: 'copy'

        input:
        file assembly from final_1
        file reads from illumina_lofreq.collect()

        output:
        file "calc_accuracy_${assembly.baseName}.tsv"

        script:
        """
        calc_accuracy.pl --ref $assembly --threads ${params.threads} $reads \
          > calc_accuracy_${assembly.baseName}.tsv
        """
    }
}

