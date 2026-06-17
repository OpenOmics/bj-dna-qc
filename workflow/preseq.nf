/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// MODULES
//

include { PUBLISH_INPUT_DATASET_WF } from '../nf-bioskryb-utils/modules/bioskryb/publish_input_dataset/main.nf'
include { CUSTOM_FASTQ_MERGE_WF } from '../nf-bioskryb-utils/modules/bioskryb/custom_fastq_merge/main.nf'
include { FastpNoQCWF; FastpQCWF } from '../nf-bioskryb-utils/modules/fastp/main.nf'
include { FASTQC } from '../nf-bioskryb-utils/modules/fastqc/main.nf'
include { SEQTK_WF } from '../nf-bioskryb-utils/modules/seqtk/sample/main.nf'
include { KRAKEN2_WF } from '../nf-bioskryb-utils/modules/kraken2/main.nf'
include { ALIGN_DEDUP_QC_SENTIEON } from '../nf-bioskryb-utils/subworkflows/align_dedup_qc/align_dedup_qc_sentieon.nf'
include { ALIGN_DEDUP_QC_STANDARD } from '../nf-bioskryb-utils/subworkflows/align_dedup_qc/align_dedup_qc_standard.nf'
include { QUALIMAP_BAMQC_WF } from '../nf-bioskryb-utils/modules/qualimap/bamqc/main.nf'
include { PRESEQ_SUBWF as PRESEQ_SUBWF} from '../nf-bioskryb-utils/modules/preseq/main.nf'
include { PRESEQ_SUBWF as PRESEQ_SUBWF_NONDEDUP} from '../nf-bioskryb-utils/modules/preseq/main.nf'
include { CUSTOM_DATA_PROCESSING_WF } from '../modules/local/custom_data_processing/main.nf'
include { CUSTOM_CALCULATE_MAPD } from '../modules/local/custom_calculate_mapd/main.nf'
include { CUSTOM_REPORT_WF } from '../modules/local/custom_report/main.nf'
include { MULTIQC_WF } from '../nf-bioskryb-utils/modules/multiqc/main.nf'
include { REPORT_VERSIONS_WF } from '../nf-bioskryb-utils/modules/bioskryb/report_tool_versions/main.nf'
include { GINKO_WF } from '../nf-bioskryb-utils/subworkflows/cnv_ginko/main.nf'
include { BAM_LORENZ_COVERAGE_WF } from '../nf-bioskryb-utils/modules/bam_lorenz_coverage/main.nf'
include { COUNT_READS_FASTQ_WF } from '../nf-bioskryb-utils/modules/bioskryb/custom_read_counts/main.nf'
include { QC_PLOTS_WF } from '../nf-bioskryb-utils/subworkflows/qc_plots/main.nf'



workflow PRESEQ_WF {
    
    take:
        ch_reads
        ch_input_csv
        ch_dummy_file
        ch_dummy_file2
        ch_mode
        ch_genome
        ch_reference
        ch_intervals
        ch_bin_size
        ch_min_ploidy
        ch_max_ploidy
        ch_min_bin_width
        ch_is_haplotype
        ch_binref
        ch_gcref
        ch_boundsref_file
        ch_seqtk_sample_seed
        ch_min_reads
        ch_n_reads
        ch_read_length
        ch_two_color_chemistry
        ch_adapter_sequence
        ch_adapter_sequence_r2
        ch_krakendb
        ch_tmp_dir
        ch_mapd_bin_size
        ch_blacklist_regions
        ch_run_sentieon
        ch_skip_seqtk
        ch_skip_kraken
        ch_skip_fastqc
        ch_skip_qualimap
        ch_skip_ginkgo
        ch_skip_mapd
        ch_publish_dir
        ch_timestamp
        ch_qcplot_config_file
        ch_enable_publish
        ch_disable_publish
    
    main:
    
        if ( ch_input_csv ) {
            
            PUBLISH_INPUT_DATASET_WF (
                                        ch_input_csv,
                                        ch_publish_dir,
                                        ch_enable_publish
                                    )
        } else {
            input_csv_dummy_file = file("input_csv_dummy_file.txt")
            ch_input_csv = Channel.fromPath( input_csv_dummy_file ).collect()
        }

        COUNT_READS_FASTQ_WF (
                                ch_reads,
                                ch_publish_dir,
                                ch_enable_publish
                            )

        COUNT_READS_FASTQ_WF.out.read_counts
        .map { sample_id, files, read_count_file -> 
            def read_count = read_count_file.text.trim().toLong()
            [sample_id, files, read_count]
        }
        .branch { read ->
            small: read[2] < ch_min_reads
            large: read[2] >= ch_min_reads
        }
        .set { branched_reads }

        ch_sample_metadata = Channel.empty()
        ch_seqtk_version = Channel.empty()
        
        if ( !ch_skip_seqtk ) {

            ch_reads_with_n_reads = branched_reads.large.map { biosampleName, reads, _read_count ->
                return [ biosampleName, reads, ch_n_reads ]
            }

            SEQTK_WF (
                            ch_reads_with_n_reads,
                            false,
                            ch_read_length,
                            ch_seqtk_sample_seed,
                            ch_publish_dir,
                            ch_disable_publish
                        )
            ch_sample_reads = SEQTK_WF.out.reads
            ch_sample_metadata = SEQTK_WF.out.metadata
            ch_seqtk_version = SEQTK_WF.out.version
        } 
        else {
            ch_sample_reads = branched_reads.large.map { biosampleName, reads, _read_count ->
                return [ biosampleName, reads ]
            }
        }

        FastpNoQCWF ( 
                ch_sample_reads, 
                ch_two_color_chemistry,
                ch_adapter_sequence,
                ch_adapter_sequence_r2,
                ch_publish_dir,
                ch_enable_publish
              )
              
        ch_fastp_report = FastpNoQCWF.out.report
        
        // FastpNoQCWF.out.reads.ifEmpty{ error "ERROR: subsample of fastq files resulted in fastq files < 10.KB in size." }
              
        ch_kraken2_report = Channel.empty()
        ch_kraken2_version = Channel.empty()

        if ( !ch_skip_kraken ) {                         
    
            KRAKEN2_WF (
                        FastpNoQCWF.out.reads,
                        ch_krakendb,
                        ch_publish_dir,
                        ch_enable_publish
                      )
                      
            ch_kraken2_report = KRAKEN2_WF.out.report
            ch_kraken2_version = KRAKEN2_WF.out.version
        }
        
        ch_fastqc_report = Channel.empty()
        ch_fastqc_version = Channel.empty()
        
        if ( !ch_skip_fastqc ) {
            FASTQC (
                        FastpNoQCWF.out.reads,
                        ch_publish_dir,
                        ch_enable_publish
                   )
            ch_fastqc_report = FASTQC.out.report
            ch_fastqc_version = FASTQC.out.version
        }

        if ( ch_run_sentieon ) {

            ALIGN_DEDUP_QC_SENTIEON (
                FastpNoQCWF.out.reads,
                ch_reference,
                ch_dummy_file,
                ch_intervals,
                ch_dummy_file2,
                ch_mode,
                ch_publish_dir,
                ch_enable_publish,
                ch_disable_publish
            )

            ch_nondedup_bam = ALIGN_DEDUP_QC_SENTIEON.out.nondedup_bam
            ch_dedup_bam = ALIGN_DEDUP_QC_SENTIEON.out.dedup_bam
            ch_dedup_metrics = ALIGN_DEDUP_QC_SENTIEON.out.dedup_metrics
            ch_metrics_with_dedup = ALIGN_DEDUP_QC_SENTIEON.out.metrics_with_dedup
            ch_metrics_with_nondedup = ALIGN_DEDUP_QC_SENTIEON.out.metrics_with_nondedup
            ch_align_dedup_qc_version = ALIGN_DEDUP_QC_SENTIEON.out.version

        } else {

            ALIGN_DEDUP_QC_STANDARD (
                FastpNoQCWF.out.reads,
                ch_reference,
                ch_dummy_file,
                ch_intervals,
                ch_dummy_file2,
                ch_mode,
                ch_publish_dir,
                ch_enable_publish,
                ch_disable_publish
            )

            ch_nondedup_bam = ALIGN_DEDUP_QC_STANDARD.out.nondedup_bam
            ch_dedup_bam = ALIGN_DEDUP_QC_STANDARD.out.dedup_bam
            ch_dedup_metrics = ALIGN_DEDUP_QC_STANDARD.out.dedup_metrics
            ch_metrics_with_dedup = ALIGN_DEDUP_QC_STANDARD.out.metrics_with_dedup
            ch_metrics_with_nondedup = ALIGN_DEDUP_QC_STANDARD.out.metrics_with_nondedup
            ch_align_dedup_qc_version = ALIGN_DEDUP_QC_STANDARD.out.version

        }
        
        ch_dedup_bam
            .collectFile( name: "bam_files.txt", newLine: true, sort: { it[0] }, storeDir: "${ch_tmp_dir}" )
                { it[0] + "\t" + "${ch_publish_dir}_${ch_timestamp}/secondary_analyses/alignment/output/" + it[1].getName() }
        
        ch_bam_lorenz_coverage_stats = Channel.empty()
        ch_bam_lorenz_coverage_version = Channel.empty()
        
        BAM_LORENZ_COVERAGE_WF (
                                ch_dedup_bam,
                                ch_publish_dir,
                                ch_enable_publish
                            )
                            
        ch_bam_lorenz_coverage_stats = BAM_LORENZ_COVERAGE_WF.out.stats
        ch_bam_lorenz_coverage_version = BAM_LORENZ_COVERAGE_WF.out.version


        ch_qualimap_report = Channel.empty()
        ch_qualimap_version = Channel.empty()

        if ( !ch_skip_qualimap ) {            
        
            QUALIMAP_BAMQC_WF ( 
                                ch_dedup_bam,
                                ch_publish_dir,
                                ch_disable_publish
                              )
            ch_qualimap_report = QUALIMAP_BAMQC_WF.out.results
            ch_qualimap_version = QUALIMAP_BAMQC_WF.out.version
        }
        
        PRESEQ_SUBWF ( 
                            ch_dedup_bam,
                            "dedup",
                            ch_publish_dir,
                            ch_disable_publish
                         )
                         
        PRESEQ_SUBWF_NONDEDUP ( 
                            ch_nondedup_bam,
                            "non_dedup",
                            ch_publish_dir,
                            ch_disable_publish
                         )
        
        combine_outputs_b = ch_dedup_bam
                                .join(ch_dedup_metrics
                                .join(ch_metrics_with_dedup
                                .join(ch_metrics_with_nondedup
                                .join(PRESEQ_SUBWF.out.coverage)
                                .join(PRESEQ_SUBWF_NONDEDUP.out.coverage)
                                                    )
                                    )
                                )
     
        
        
        // CNV CALLING = GINKGO
        ch_ginkgo_version = Channel.empty()
        ch_ginkgo_stats = Channel.empty()
        ch_bedtools_version = Channel.empty()
        if( !ch_skip_ginkgo && (ch_genome == "GRCh38" || ch_genome == "GRCm39" || ch_genome == "ARSUCD2") ) {
            GINKO_WF(
                    ch_dedup_bam,
                    ch_bin_size,
                    ch_binref,
                    ch_gcref,
                    ch_boundsref_file,
                    ch_min_ploidy,
                    ch_max_ploidy,
                    ch_min_bin_width,
                    ch_is_haplotype,
                    ch_publish_dir,
                    ch_enable_publish,
                    ch_disable_publish
            )

            ch_ginkgo_version = GINKO_WF.out.ginkgo_version
            ch_ginkgo_stats = GINKO_WF.out.ginkgo_stats
            ch_bedtools_version = GINKO_WF.out.bedtools_version

        }
        
        CUSTOM_DATA_PROCESSING_WF ( 
                                  combine_outputs_b,
                                  ch_reference,
                                  ch_publish_dir,
                                  ch_disable_publish
                                )

        ch_custom_calculate_mapd_metrics = channel.empty()
        ch_custom_calculate_mapd_version = channel.empty()

        if ( !ch_skip_mapd ) {
            
            CUSTOM_CALCULATE_MAPD ( 
                                        ch_dedup_bam,
                                        ch_mapd_bin_size,
                                        ch_blacklist_regions,
                                        ch_reference,
                                        ch_publish_dir,
                                        ch_disable_publish
                                     )
            
            ch_custom_calculate_mapd_metrics = CUSTOM_CALCULATE_MAPD.out.metrics
            ch_custom_calculate_mapd_version = CUSTOM_CALCULATE_MAPD.out.version
        }
        
        
        CUSTOM_REPORT_WF ( 
                            CUSTOM_DATA_PROCESSING_WF.out.metrics
                            .collect()
                            .combine( ch_sample_metadata.collect().ifEmpty([]) )
                            .combine( ch_fastp_report.collect().ifEmpty([]) )
                            .combine( ch_bam_lorenz_coverage_stats.collect().ifEmpty([]) )
                            .combine( ch_custom_calculate_mapd_metrics.collect().ifEmpty([]) )
                            .combine( ch_ginkgo_stats.ifEmpty([]) ),
                            COUNT_READS_FASTQ_WF.out.combined_read_counts,
                            ch_min_reads,
                            ch_publish_dir,
                            ch_enable_publish
                        )

        ch_tool_versions = ch_seqtk_version.take(1).ifEmpty([])
                            .combine(FastpNoQCWF.out.version.take(1))
                            .combine(ch_fastqc_version.take(1).ifEmpty([]))
                            .combine(ch_align_dedup_qc_version.collect())
                            .combine(ch_qualimap_version.take(1).ifEmpty([]) )
                            .combine(PRESEQ_SUBWF.out.preseqBam2mr_version.take(1))
                            .combine(PRESEQ_SUBWF.out.preseqExtrap_version.take(1))
                            .combine(ch_kraken2_version.take(1).ifEmpty([]) )
                            .combine(ch_bam_lorenz_coverage_version.take(1).ifEmpty([]))
                            .combine((ch_ginkgo_version ?: Channel.empty()).ifEmpty([]))
                            .combine(ch_bedtools_version.take(1).ifEmpty([]))
                                                    
        REPORT_VERSIONS_WF(
                            ch_tool_versions,
                            ch_publish_dir,
                            ch_enable_publish
                          )
        qc_plots_composition_jpg = channel.empty()
        qc_plots_cnv_quadrants_jpg = channel.empty()
        ch_custom_report_cnv = channel.empty()
        
        // Run the QC_Plots for all Ginkgo runs
        if( !ch_skip_ginkgo && (ch_genome == "GRCh38" || ch_genome == "GRCm39" || ch_genome == "ARSUCD2") ) {
            QC_PLOTS_WF(
                GINKO_WF.out.RDS,
                GINKO_WF.out.SEGCOPY,
                CUSTOM_REPORT_WF.out.all_metrics,
                CUSTOM_REPORT_WF.out.selected_metrics,
                ch_input_csv,
                ch_qcplot_config_file,
                ch_publish_dir,
                ch_enable_publish
            )
            
            qc_plots_composition_jpg = QC_PLOTS_WF.out.composition_jpg.ifEmpty([])
            qc_plots_cnv_quadrants_jpg = QC_PLOTS_WF.out.cnv_quadrants_jpg.ifEmpty([])
            ch_custom_report_cnv = QC_PLOTS_WF.out.allmetrics_with_cnv.mix( QC_PLOTS_WF.out.selectedmetrics_with_cnv ).collect()
            ch_fallback_report = CUSTOM_REPORT_WF.out.mqc
            // Combine both and use the first available
            ch_custom_report = ch_custom_report_cnv.concat(ch_fallback_report).first() // Use CNV-enhanced report if available, otherwise fall back to original custom report
        }
        else {
            ch_custom_report = CUSTOM_REPORT_WF.out.mqc
        }

        // Only include QC plots outputs if they exist (i.e., QC_PLOTS_WF succeeded)
        collect_mqc = ch_custom_report
                        .combine( ch_fastp_report.collect() )
                        .combine( ch_fastqc_report.collect().ifEmpty([]) )
                        .combine( ch_kraken2_report.collect().ifEmpty([]) )
                        .combine( ch_qualimap_report.collect().ifEmpty([]) )
                        .combine( REPORT_VERSIONS_WF.out.versions.collect().ifEmpty([]))
                        .combine(qc_plots_composition_jpg.collect().ifEmpty([]))
                        .combine(qc_plots_cnv_quadrants_jpg.collect().ifEmpty([]))


    emit:
        dedup_bam = ch_dedup_bam
        multiqc_input = collect_mqc
        // params_meta = params_meta
  
}


// OnComplete
workflow.onComplete {
    
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]
    output["output"]["bam"]         = [:]

    bam_outfile = file("$params.tmp_dir/bam_files.txt")
    bam_outfile_lines = bam_outfile.readLines()
    for ( bam_line : bam_outfile_lines ) {
        def (sample_name, bam_path) = bam_line.split('\t')
        output["output"]["bam"][sample_name] = [:]
        output["output"]["bam"][sample_name]["bam"] = bam_path
    }
    
    def output_json = groovy.json.JsonOutput.toJson(output)
    def output_json_pretty = groovy.json.JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    println(output_json_pretty)

    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError {
    
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]
    output["output"]["bam"]         = [:]
    output["output"]["vcf"]         = [:]
    
    def output_json = groovy.json.JsonOutput.toJson(output)
    def output_json_pretty = groovy.json.JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    println(output_json_pretty)

    def subject = """\
        [nf-preseq-pipeline] FAILED: ${workflow.runName}
        """

    def msg = """\

        Pipeline execution summary 
        --------------------------------
        Script name       : ${workflow.scriptName ?: '-'}
        Script ID         : ${workflow.scriptId ?: '-'}
        Workflow session  : ${workflow.sessionId}
        Workflow repo     : ${workflow.repository ?: '-' }
        Workflow revision : ${workflow.repository ? "$workflow.revision ($workflow.commitId)" : '-'}
        Workflow profile  : ${workflow.profile ?: '-'}
        Workflow cmdline  : ${workflow.commandLine ?: '-'}
        Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
        Error Report      : ${workflow.errorReport}
        """
        .stripIndent()
    
    log.info ( msg )
    
    if ( "${params.email_on_fail}" && workflow.exitStatus != 0 ) {
        sendMail(to: "${params.email_on_fail}", subject: subject, body: msg)
    }

    println( "\nPipeline failed.\n\n" )
}
