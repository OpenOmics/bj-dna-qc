nextflow.enable.dsl=2
params.timestamp = ""

process PICARD_METRICS {
    tag "${sample_name}"
    publishDir "${params.publish_dir}_${params.timestamp}/secondary_analyses/metrics/${sample_name}/alignment_stats", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(bam), path(bai), path(recal_table_file)
    path fasta_ref
    path base_metrics_intervals
    path wgs_or_target_intervals
    val mode
    val type
    val(publish_dir)
    val(enable_publish)


    output:
    tuple val(sample_name), file("*metrics*"), emit: metrics_tuple
    path "*metrics*", emit: metrics
    tuple val(sample_name), file("*.${type}.alignmentstat_metrics.txt"), emit: alignment_metrics
    path("picard_version.yml"), emit: version
    
    script:
    def bqsr = recal_table_file.name == "dummy_file.txt" ? "" : "--bqsr-recal-file ${recal_table_file}"
    if (mode == 'exome') {
        """
        samtools view -b -L ${wgs_or_target_intervals} ${bam} -o ${sample_name}.${type}.interval.bam

        picard CollectMultipleMetrics \
            -I ${sample_name}.${type}.interval.bam \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type} \
            --PROGRAM null \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM MeanQualityByCycle
        
        ln -s ${sample_name}.${type}.gc_bias.summary_metrics ${sample_name}.${type}.gcbias_summary.metrics.txt
        ln -s ${sample_name}.${type}.gc_bias.detail_metrics ${sample_name}.${type}.gcbias.metrics.txt
        ln -s ${sample_name}.${type}.alignment_summary_metrics ${sample_name}.${type}.alignmentstat_metrics.txt
        ln -s ${sample_name}.${type}.insert_size_metrics ${sample_name}.${type}.insertsizemetricalgo.metrics.txt
        ln -s ${sample_name}.${type}.quality_by_cycle_metrics ${sample_name}.${type}.meanqualitybycycle.metrics.txt

        samtools bedcov -j -H ${wgs_or_target_intervals} ${bam} > ${sample_name}.${type}.cov_metrics.sample_interval_summary

        picard BedToIntervalList \
            -I ${wgs_or_target_intervals} \
            -O ${wgs_or_target_intervals}.interval_list \
            -SD ${fasta_ref}/genome.fa

        #picard CollectHsMetrics \
        #    -I ${bam} \
        #    -R ${fasta_ref}/genome.fa \
        #    -O ${sample_name}.${type}.hsmetricalgo.metrics.txt \
        #    --TARGET_INTERVALS ${wgs_or_target_intervals}.interval_list \
        #    --BAIT_INTERVALS ${wgs_or_target_intervals}.interval_list

        export PICARD_VER=\$(echo \$(picard CollectMultipleMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        echo Picard: \$PICARD_VER > picard_version.yml
        """
    } else if (mode == 'wgs') {
        """
        echo "${bqsr}"
        echo "${recal_table_file.name}"
        
        samtools view -b -L ${base_metrics_intervals} ${bam} -o ${sample_name}.${type}.interval.bam

        picard CollectMultipleMetrics \
            -I ${sample_name}.${type}.interval.bam \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type} \
            --PROGRAM null \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM MeanQualityByCycle
        
        ln -s ${sample_name}.${type}.gc_bias.summary_metrics ${sample_name}.${type}.gcbias_summary.metrics.txt
        ln -s ${sample_name}.${type}.gc_bias.detail_metrics ${sample_name}.${type}.gcbias.metrics.txt
        ln -s ${sample_name}.${type}.alignment_summary_metrics ${sample_name}.${type}.alignmentstat_metrics.txt
        ln -s ${sample_name}.${type}.insert_size_metrics ${sample_name}.${type}.insertsizemetricalgo.metrics.txt
        ln -s ${sample_name}.${type}.quality_by_cycle_metrics ${sample_name}.${type}.meanqualitybycycle.metrics.txt

        samtools bedcov -j -H ${base_metrics_intervals} ${bam} > ${sample_name}.${type}.cov_metrics.sample_interval_summary

        picard BedToIntervalList \
            -I ${wgs_or_target_intervals} \
            -O ${wgs_or_target_intervals}.interval_list \
            -SD ${fasta_ref}/genome.fa

        picard CollectWgsMetrics \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type}.wgsmetricsalgo.metrics.txt \
            --INTERVALS ${wgs_or_target_intervals}.interval_list
        
        export PICARD_VER=\$(echo \$(picard CollectMultipleMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        echo Picard: \$PICARD_VER > picard_version.yml
        """
    } else {
        
        """
        samtools view -b -L ${base_metrics_intervals} ${bam} -o ${sample_name}.${type}.interval.bam

        picard CollectMultipleMetrics \
            -I ${sample_name}.${type}.interval.bam \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type} \
            --PROGRAM null \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM MeanQualityByCycle

        ln -s ${sample_name}.${type}.gc_bias.summary_metrics ${sample_name}.${type}.gcbias_summary.metrics.txt
        ln -s ${sample_name}.${type}.gc_bias.detail_metrics ${sample_name}.${type}.gcbias.metrics.txt
        ln -s ${sample_name}.${type}.alignment_summary_metrics ${sample_name}.${type}.alignmentstat_metrics.txt
        ln -s ${sample_name}.${type}.insert_size_metrics ${sample_name}.${type}.insertsizemetricalgo.metrics.txt
        ln -s ${sample_name}.${type}.quality_by_cycle_metrics ${sample_name}.${type}.meanqualitybycycle.metrics.txt

        samtools bedcov -j -H ${base_metrics_intervals} ${bam} > ${sample_name}.${type}.cov_metrics.sample_interval_summary

        picard BedToIntervalList \
            -I ${base_metrics_intervals} \
            -O ${base_metrics_intervals}.interval_list \
            -SD ${fasta_ref}/genome.fa
        
        picard CollectWgsMetrics \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type}.wgsmetricsalgo.metrics.txt \
            --INTERVALS ${base_metrics_intervals}.interval_list

        export PICARD_VER=\$(echo \$(picard CollectMultipleMetrics --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
        echo Picard: \$PICARD_VER > picard_version.yml
        """
    }
}



workflow PICARD_METRICS_WF {
    
    take:
        ch_bam
        ch_reference
        ch_base_metrics_intervals
        ch_wgs_or_target_intervals
        ch_mode
        ch_type
        ch_publish_dir
        ch_enable_publish
        
    main:
        PICARD_METRICS ( 
            ch_bam,
            ch_reference,
            ch_base_metrics_intervals,
            ch_wgs_or_target_intervals,
            ch_mode,
            ch_type,
            ch_publish_dir,
            ch_enable_publish
        )
            
                              
    emit:
        metrics_tuple = PICARD_METRICS.out.metrics_tuple
        metrics = PICARD_METRICS.out.metrics
        alignment_metrics = PICARD_METRICS.out.alignment_metrics
        version = PICARD_METRICS.out.version
}

include { CUSTOM_METRICS_MERGE } from '../../bioskryb/custom_metrics_merge/main.nf'


workflow {
    
    ch_dummy_file = Channel.fromPath(params.dummy_file, checkIfExists: true).collect()
    
    if (params.bam != "") {
        ch_bam_raw = Channel.fromFilePairs(params.bam, size: -1)
        ch_bam_raw.map{ it -> it.flatten().collect() }
        ch_bam = ch_bam_raw.combine(ch_dummy_file)
    } else if(params.input_csv != "") {
        ch_bam_raw = Channel.fromPath(params.input_csv).splitCsv(header:true)
                                .map { row -> [ row.biosampleName, row.bam, row.bam + ".bai" ] }
        ch_bam = ch_bam_raw.combine(ch_dummy_file)
    }
    
    ch_bam.view()
    ch_bam.ifEmpty{ exit 1, "ERROR: No BAM files specified either via --bam or --input_csv" }
    
    PICARD_METRICS_WF( 
        ch_bam,
        params.reference,
        params.base_metrics_intervals,
        params.wgs_or_target_intervals,
        params.mode,
        "dedup",
        params.publish_dir,
        params.enable_publish
    )
    
    CUSTOM_METRICS_MERGE(
        PICARD_METRICS_WF.out.metrics.collect(),
        params.publish_dir,
        params.enable_publish
    )
}
