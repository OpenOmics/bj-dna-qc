nextflow.enable.dsl=2

include { BWA_WF } from '../../modules/bwa/mem/main.nf'
include { PICARD_MARKDUPLICATES_WF } from '../../modules/picard/markduplicates/main.nf'
include { GATK4_METRICS_WF as GATK4_METRICS_WITH_DEDUP_WF } from '../../modules/gatk4/metrics/main.nf'
include { GATK4_METRICS_WF as GATK4_METRICS_WITH_NONDEDUP_WF } from '../../modules/gatk4/metrics/main.nf'


workflow ALIGN_DEDUP_QC_STANDARD {
    take:
        ch_reads
        ch_reference
        ch_dummy_file
        ch_base_metrics_intervals
        ch_wgs_or_target_intervals
        ch_mode
        ch_publish_dir
        ch_enable_publish
        ch_disable_publish

    main:
        BWA_WF (
            ch_reads,
            ch_reference,
            ch_publish_dir,
            ch_disable_publish
        )

        PICARD_MARKDUPLICATES_WF (
            BWA_WF.out.bam,
            ch_reference,
            ch_publish_dir,
            ch_enable_publish
        )

        GATK4_METRICS_WITH_DEDUP_WF (
            PICARD_MARKDUPLICATES_WF.out.bam.combine(ch_dummy_file),
            ch_reference,
            ch_base_metrics_intervals,
            ch_wgs_or_target_intervals,
            ch_mode,
            "dedup",
            ch_publish_dir,
            ch_enable_publish
        )

        GATK4_METRICS_WITH_NONDEDUP_WF (
            BWA_WF.out.bam.combine(ch_dummy_file),
            ch_reference,
            ch_base_metrics_intervals,
            ch_wgs_or_target_intervals,
            ch_mode,
            "nondedup",
            ch_publish_dir,
            ch_enable_publish
        )

        ch_versions = Channel.empty()
        ch_versions = ch_versions.mix(BWA_WF.out.version.first())
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES_WF.out.version.first())
        ch_versions = ch_versions.mix(GATK4_METRICS_WITH_DEDUP_WF.out.version.first())

    emit:
        nondedup_bam = BWA_WF.out.bam
        dedup_bam = PICARD_MARKDUPLICATES_WF.out.bam
        dedup_metrics = PICARD_MARKDUPLICATES_WF.out.metrics
        metrics_with_dedup = GATK4_METRICS_WITH_DEDUP_WF.out.metrics_tuple
        metrics_with_nondedup = GATK4_METRICS_WITH_NONDEDUP_WF.out.metrics_tuple
        version = ch_versions
}
