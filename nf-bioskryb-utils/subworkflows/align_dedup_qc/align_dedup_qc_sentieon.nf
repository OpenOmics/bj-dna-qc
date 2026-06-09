nextflow.enable.dsl=2

include { SENTIEON_BWA_WF } from '../../modules/sentieon/bwa/mem/main.nf'
include { SENTIEON_DRIVER_LOCUSCOLLECTOR_WF } from '../../modules/sentieon/driver/locuscollector/main.nf'
include { SENTIEON_DRIVER_DEDUP_WF } from '../../modules/sentieon/driver/dedup/main.nf'
include { SENTIEON_DRIVER_METRICS_WF as SENTIEON_DRIVER_METRICS_WITH_DEDUP_WF } from '../../modules/sentieon/driver/metrics/main.nf'
include { SENTIEON_DRIVER_METRICS_WF as SENTIEON_DRIVER_METRICS_WITH_NONDEDUP_WF } from '../../modules/sentieon/driver/metrics/main.nf'


workflow ALIGN_DEDUP_QC_SENTIEON {
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
        SENTIEON_BWA_WF (
            ch_reads,
            ch_reference,
            ch_publish_dir,
            ch_disable_publish
        )

        SENTIEON_DRIVER_LOCUSCOLLECTOR_WF (
            SENTIEON_BWA_WF.out.bam,
            ch_reference,
            ch_publish_dir,
            ch_disable_publish
        )

        combine_outputs_a = SENTIEON_BWA_WF.out.bam.join(SENTIEON_DRIVER_LOCUSCOLLECTOR_WF.out.locuscollector_score)

        SENTIEON_DRIVER_DEDUP_WF (
            combine_outputs_a,
            ch_reference,
            ch_publish_dir,
            ch_enable_publish
        )

        SENTIEON_DRIVER_METRICS_WITH_DEDUP_WF (
            SENTIEON_DRIVER_DEDUP_WF.out.bam.combine(ch_dummy_file),
            ch_reference,
            ch_base_metrics_intervals,
            ch_wgs_or_target_intervals,
            ch_mode,
            "dedup",
            ch_publish_dir,
            ch_enable_publish
        )

        SENTIEON_DRIVER_METRICS_WITH_NONDEDUP_WF (
            SENTIEON_BWA_WF.out.bam.combine(ch_dummy_file),
            ch_reference,
            ch_base_metrics_intervals,
            ch_wgs_or_target_intervals,
            ch_mode,
            "nondedup",
            ch_publish_dir,
            ch_enable_publish
        )

        ch_versions = Channel.empty()
        ch_versions = ch_versions.mix(SENTIEON_BWA_WF.out.version.first())
        ch_versions = ch_versions.mix(SENTIEON_DRIVER_LOCUSCOLLECTOR_WF.out.version.first())
        ch_versions = ch_versions.mix(SENTIEON_DRIVER_DEDUP_WF.out.version.first())
        ch_versions = ch_versions.mix(SENTIEON_DRIVER_METRICS_WITH_DEDUP_WF.out.version.first())

    emit:
        nondedup_bam = SENTIEON_BWA_WF.out.bam
        dedup_bam = SENTIEON_DRIVER_DEDUP_WF.out.bam
        dedup_metrics = SENTIEON_DRIVER_DEDUP_WF.out.metrics
        metrics_with_dedup = SENTIEON_DRIVER_METRICS_WITH_DEDUP_WF.out.metrics_tuple
        metrics_with_nondedup = SENTIEON_DRIVER_METRICS_WITH_NONDEDUP_WF.out.metrics_tuple
        version = ch_versions
}
