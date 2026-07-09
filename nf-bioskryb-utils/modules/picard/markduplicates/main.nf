nextflow.enable.dsl=2
params.timestamp = ""

process PICARD_MARKDUPLICATES {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/secondary_analyses/alignment", enabled:"$enable_publish", pattern: "output/${sample_name}.bam*"

    input:
    tuple val(sample_name), path(bams), path(bais)
    path fasta_ref
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), path("output/${sample_name}.bam"), path("output/${sample_name}.bam.bai"), emit: bam
    tuple val(sample_name), path("*.dedup_metrics.txt"), emit: metrics
    path("*.dedup_metrics.txt"), emit: mqc_metrics
    path("output/${sample_name}.bam"), emit: bam_to_compress
    path("picard_markduplicates_version.yml"), emit: version

    script:

    """
    mkdir output;

    picard MarkDuplicates \
        --INPUT ${bams.join(" --INPUT ")} \
        --OUTPUT output/${sample_name}.bam \
        --METRICS_FILE ${sample_name}.dedup_metrics.txt \
        --REFERENCE_SEQUENCE ${fasta_ref}/genome.fa \
        --REMOVE_DUPLICATES \
        --CREATE_INDEX
    
    mv output/${sample_name}.bai output/${sample_name}.bam.bai

    export PICARD_VER=\$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    echo Picard: \$PICARD_VER > picard_markduplicates_version.yml
    """
}

workflow PICARD_MARKDUPLICATES_WF{
    take:
        ch_bam
        ch_reference
        ch_publish_dir
        ch_enable_publish
        
    main:
        PICARD_MARKDUPLICATES ( 
                                ch_bam,
                                ch_reference,
                                ch_publish_dir,
                                ch_enable_publish
                              )
    emit:
        bam = PICARD_MARKDUPLICATES.out.bam
        metrics = PICARD_MARKDUPLICATES.out.metrics
        mqc_metrics = PICARD_MARKDUPLICATES.out.mqc_metrics
        bam_to_compress = PICARD_MARKDUPLICATES.out.bam_to_compress
        version = PICARD_MARKDUPLICATES.out.version
    
}

workflow{
    
    ch_bam = Channel.fromFilePairs(params.bam_dir + "/*{.bam,.bam.bai}")
    
    PICARD_MARKDUPLICATES_WF ( 
                                ch_bam,
                                params.reference,
                                params.publish_dir,
                                params.enable_publish
                             )
    
}
