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
    tuple val(sample_name), file("*metrics.txt*"), emit: metrics_tuple
    path "*metrics.txt*", emit: metrics
    tuple val(sample_name), file("*.${type}.alignmentstat_metrics.txt"), emit: alignment_metrics
    path("gatk4_qcmetrics_version.yml"), emit: version
    
    script:
    def bqsr = recal_table_file.name == "dummy_file.txt" ? "" : "--bqsr-recal-file ${recal_table_file}"
    if (mode == 'exome') {
        """
        gatk CollectMultipleMetrics \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type} \
            -L ${wgs_or_target_intervals} \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM MeanQualityByCycle

        ln -s ${sample_name}.${type}.gc_bias.summary_metrics ${sample_name}.${type}.gcbias_summary.metrics.txt
        ln -s ${sample_name}.${type}.gc_bias.detail_metrics ${sample_name}.${type}.gcbias.metrics.txt
        ln -s ${sample_name}.${type}.alignment_summary_metrics ${sample_name}.${type}.alignmentstat_metrics.txt
        ln -s ${sample_name}.${type}.insert_size_metrics ${sample_name}.${type}.insertsizemetricalgo.metrics.txt
        ln -s ${sample_name}.${type}.quality_by_cycle_metrics ${sample_name}.${type}.meanqualitybycycle.metrics.txt

        gatk DepthOfCoverage \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -L ${wgs_or_target_intervals} \
            -O ${sample_name}.${type}.cov_metrics

        #gatk CollectHsMetrics \
        #    -I ${bam} \
        #    -R ${fasta_ref}/genome.fa \
        #    -O ${sample_name}.${type}.hsmetricalgo.metrics.txt \
        #    --TARGET_INTERVALS ${wgs_or_target_intervals} \
        #    --BAIT_INTERVALS ${wgs_or_target_intervals}

        # These files are not created in exome. creating dummy files?
        touch ${sample_name}.wgsmetricsalgo.metrics.txt

        export GATK4_VER=\$(echo \$(gatk --version 2>&1) | sed -e 's/.*(GATK) //; s/Version: //g')
        echo GATK4: \$GATK4_VER > gatk4_qcmetrics_version.yml
        """
    } else if (mode == 'wgs') {
        """
        echo "${bqsr}"
        echo "${recal_table_file.name}"
        
        gatk CollectMultipleMetrics \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type} \
            --INTERVALS ${base_metrics_intervals} \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM MeanQualityByCycle
        
        ln -s ${sample_name}.${type}.gc_bias.summary_metrics ${sample_name}.${type}.gcbias_summary.metrics.txt
        ln -s ${sample_name}.${type}.gc_bias.detail_metrics ${sample_name}.${type}.gcbias.metrics.txt
        ln -s ${sample_name}.${type}.alignment_summary_metrics ${sample_name}.${type}.alignmentstat_metrics.txt
        ln -s ${sample_name}.${type}.insert_size_metrics ${sample_name}.${type}.insertsizemetricalgo.metrics.txt
        ln -s ${sample_name}.${type}.quality_by_cycle_metrics ${sample_name}.${type}.meanqualitybycycle.metrics.txt

        gatk DepthOfCoverage \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -L ${base_metrics_intervals} \
            -O ${sample_name}.${type}.cov_metrics \
            --omit-depth-output-at-each-base \
            --omit-locus-table \
            --omit-per-sample-statistics

        gatk CollectWgsMetrics \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type}.wgsmetricsalgo.metrics.txt \
            --INTERVALS ${wgs_or_target_intervals}
                
        # FROM CROMWELL - these files are not created in wgs. creating dummy files?
        # touch ${sample_name}.${type}.coveragemetrics.metrics.sample_summary
        touch ${sample_name}.${type}.hsmetricalgo.metrics.txt
        
        export GATK4_VER=\$(echo \$(gatk --version 2>&1) | sed -e 's/.*(GATK) //; s/Version: //g')
        echo GATK4: \$GATK4_VER > gatk4_qcmetrics_version.yml
        """
    } else {
        
        """
        gatk BedToIntervalList \
            -I ${base_metrics_intervals} \
            -O ${base_metrics_intervals}.interval_list \
            -SD ${fasta_ref}/genome.fa

        gatk CollectMultipleMetrics \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type} \
            --INTERVALS ${base_metrics_intervals}.interval_list \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM MeanQualityByCycle

        ln -s ${sample_name}.${type}.gc_bias.summary_metrics ${sample_name}.${type}.gcbias_summary.metrics.txt
        ln -s ${sample_name}.${type}.gc_bias.detail_metrics ${sample_name}.${type}.gcbias.metrics.txt
        ln -s ${sample_name}.${type}.alignment_summary_metrics ${sample_name}.${type}.alignmentstat_metrics.txt
        ln -s ${sample_name}.${type}.insert_size_metrics ${sample_name}.${type}.insertsizemetricalgo.metrics.txt
        ln -s ${sample_name}.${type}.quality_by_cycle_metrics ${sample_name}.${type}.meanqualitybycycle.metrics.txt

        gatk DepthOfCoverage \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -L ${base_metrics_intervals} \
            -O ${sample_name}.${type}.cov_metrics \
            --omit-depth-output-at-each-base \
            --omit-locus-table \
            --omit-per-sample-statistics
        
        gatk CollectWgsMetrics \
            -I ${bam} \
            -R ${fasta_ref}/genome.fa \
            -O ${sample_name}.${type}.wgsmetricsalgo.metrics.txt \
            --INTERVALS ${base_metrics_intervals}.interval_list

        export GATK4_VER=\$(echo \$(gatk --version 2>&1) | sed -e 's/.*(GATK) //; s/Version: //g')
        echo GATK4: \$GATK4_VER > gatk4_qcmetrics_version.yml
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
