params.publish_dir = ""
params.timestamp = ""

process CUSTOM_REPORT {
    tag 'report'
    publishDir "${params.publish_dir}_${params.timestamp}/secondary_analyses/metrics", enabled:"$enable_publish"

    
    input:
    path(files)
    path(read_counts_csv)
    val(min_reads)
    val(publish_dir)
    val(enable_publish)


    output:
    path "nf-preseq-pipeline*", emit: mqc
    path "nf-preseq-pipeline_all_metrics*", emit: all_metrics
    path "nf-preseq-pipeline_selected_metrics*", emit: selected_metrics
    path("custom_report_version.yml"), emit: version

    
    script:
    """
    python3 /scripts/custom_report.py -s '*metrics.txt' -m '*MAPD' -q '*no_qc_fastp.json' -t '*sample_metadata.json' -l '*lorenzstats.tsv' -i '$read_counts_csv' -r '$min_reads' -o 'Summary'


    echo custom_report: v0.0.1 > custom_report_version.yml
    
    """
    
}

workflow CUSTOM_REPORT_WF{
    take:
        ch_metrics
        ch_read_counts_csv
        ch_min_reads
        ch_publish_dir
        ch_enable_publish
    main:
        CUSTOM_REPORT ( 
                        ch_metrics,
                        ch_read_counts_csv,
                        ch_min_reads,
                        ch_publish_dir,
                        ch_enable_publish
                      )
    emit:
        mqc = CUSTOM_REPORT.out.mqc
        all_metrics = CUSTOM_REPORT.out.all_metrics
        selected_metrics = CUSTOM_REPORT.out.selected_metrics
        version = CUSTOM_REPORT.out.version
}

workflow{
    
    ch_files = Channel.fromPath(params.metrics_file)
    ch_read_counts_csv = Channel.fromPath("NO_FILE")
    ch_min_reads = 1000
    
    CUSTOM_REPORT_WF ( ch_files, ch_read_counts_csv, ch_min_reads, params.publish_dir, params.enable_publish )
    
}
