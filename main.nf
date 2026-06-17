include { printHeader ; helpMessage } from './help'
include { PRESEQ_WF   } from './workflow/preseq.nf'
include { PSEUDO_BULK_WF } from './nf-bioskryb-utils/subworkflows/pseudobulk_sc_wf/main.nf'
include { SENTIEON_DNASCOPE } from './nf-bioskryb-utils/modules/sentieon/driver/dnascope/main.nf'
include { BCFTOOLS_ISEC } from './nf-bioskryb-utils/modules/bcftools/filter_isec/main.nf'
include { SIGPROFILERGENERATEMATRIX_WF } from './nf-bioskryb-utils/modules/sigprofilermatrixgenerator/main.nf'
include { GOOGLE_DEEPVARIANT_WF } from './nf-bioskryb-utils/subworkflows/google_deepvariant_wf/main.nf'
include { MULTIQC_WF } from './nf-bioskryb-utils/modules/multiqc/main.nf'


workflow {

    if (params.help) {
        helpMessage()
        System.exit(0)
    }

    printHeader()

    // Check if publishDir specified

    if (params.publish_dir == "") {
        error "ERROR: publish_dir is not defined.\nPlease add --publish_dir </path/to/output_dir> to specify where the pipeline outputs will be stored."
    }

    if (params.reads) {

        ch_reads = Channel.fromFilePairs(params.reads, size: -1, checkExists: true)
        ch_reads.ifEmpty { error "ERROR: cannot find any fastq files matching the pattern: ${params.reads}\nMake sure that the input file exists!" }
    }
    else if (params.input_csv) {

        ch_reads = Channel
            .fromPath(params.input_csv)
            .splitCsv(header: true)
            .map { row -> [row.biosampleName, [row.read1, row.read2]] }
        ch_reads.ifEmpty { error "ERROR: Input csv file is empty." }
    }

    ch_reads.view()
    ch_dummy_file = Channel.fromPath("${projectDir}/assets/dummy_file.txt", checkIfExists: true).collect()
    ch_dummy_file2 = Channel.fromPath("${projectDir}/assets/dummy_file2.txt", checkIfExists: true).collect()
    ch_binref = Channel.fromPath(params.ginko_ref_dir + "variable_" + params.bin_size + "_" + params.read_length + "_bwa")
    ch_gcref = Channel.fromPath(params.ginko_ref_dir + "GC_variable_" + params.bin_size + "_" + params.read_length + "_bwa")
    ch_boundsref_file = Channel.fromPath(params.ginko_ref_dir + "bounds_variable_" + params.bin_size + "_" + params.read_length + "_bwa")

    if ( params.genome == 'GRCh38' ) {
            dnascope_model  = params.genomes [ params.genome ] [ params.platform ] [ params.dnascope_model_selection ] [ 'wgs' ] [ 'dnascope_model' ]
            vcfeval_baseline_vcf         = params.genomes[params.genome] [ params.giab_reference_name ] [ 'vcfeval_baseline_vcf' ]
            vcfeval_baseline_vcf_index   = params.genomes[params.genome] [ params.giab_reference_name ] [ 'vcfeval_baseline_vcf_index' ]
            sigprofilermatrixgenerator_reference  = params.genomes[params.genome]['sigprofilermatrixgenerator_reference']
            sigprofilermatrixgenerator_interval   = params.genomes [ params.genome ] [ 'sigprofilermatrixgenerator_interval' ]
    }

    if (!params.is_bam){
        PRESEQ_WF(
            ch_reads,
            params.input_csv,
            ch_dummy_file,
            ch_dummy_file2,
            params.mode,
            params.genome,
            params.reference,
            params.intervals,
            params.bin_size,
            params.min_ploidy,
            params.max_ploidy,
            params.min_bin_width,
            params.is_haplotype,
            ch_binref,
            ch_gcref,
            ch_boundsref_file,
            params.seqtk_sample_seed,
            params.min_reads,
            params.n_reads,
            params.read_length,
            params.two_color_chemistry,
            params.adapter_sequence,
            params.adapter_sequence_r2,
            params.krakendb,
            params.tmp_dir,
            params.mapd_bin_size,
            params.blacklist_regions,
            params.run_sentieon,
            params.skip_seqtk,
            params.skip_kraken,
            params.skip_fastqc,
            params.skip_qualimap,
            params.skip_ginkgo,
            params.skip_mapd,
            params.publish_dir,
            params.timestamp,
            params.plot_qc_config,
            params.enable_publish,
            params.disable_publish,
        )
    }

    // Mutational Signatures Profiling
    if (!params.skip_sigprofile && params.genome == 'GRCh38') {
        if (!params.is_bam){
            ch_reads_with_groups = Channel
                .fromPath(params.input_csv)
                .splitCsv(header: true)
                .map { row -> [row.biosampleName, [row.read1, row.read2], row.groups] }
                .branch {
                    with_groups: it[2] != null && it[2] != "" && it[2] != "None"
                    without_groups: true
                }

            def samples_with_groups = ch_reads_with_groups.with_groups

            pseudobulk_input = samples_with_groups
                .map { sample_name, _fastqs, groups -> [sample_name, groups] }
                .join(PRESEQ_WF.out.dedup_bam)
                .map { sample_name, groups, bam, bai -> 
                    [sample_name, bam, bai, groups, params.n_reads]
                }
        } else {
            pseudobulk_input = Channel
                .fromPath(params.input_csv)
                .splitCsv(header: true)
                .map { row -> [row.biosampleName, row.bam, row.bam + '.bai', row.groups, params.n_reads] }
        }

        PSEUDO_BULK_WF (
                            pseudobulk_input,
                            params.genome,
                            params.reference,
                            params.dbsnp,
                            params.dbsnp_index,
                            params.mills,
                            params.mills_index,
                            params.onekg_omni,
                            params.onekg_omni_index,
                            params.publish_dir,
                            params.enable_publish
                        )

        SENTIEON_DNASCOPE (
                            PSEUDO_BULK_WF.out.pseudo_bam,
                            params.reference,
                            params.calling_intervals_filename,
                            params.dbsnp,
                            params.dbsnp_index,
                            dnascope_model,
                            params.pcrfree,
                            params.ploidy,
                            "gvcf",
                            params.publish_dir,
                            params.enable_publish
                        )

        ch_vcf = SENTIEON_DNASCOPE.out.vcf.map
            { sample_name, vcf -> [sample_name, vcf[0], vcf[1]] }

        
        BCFTOOLS_ISEC (
                        ch_vcf,
                        vcfeval_baseline_vcf,
                        vcfeval_baseline_vcf_index,
                        params.reference,
                        params.dbsnp,
                        params.dbsnp_index,
                        params.publish_dir,
                        params.enable_publish
                    )

        SIGPROFILERGENERATEMATRIX_WF (    
            BCFTOOLS_ISEC.out.filtered_vcf.map { sample_name, vcf -> [sample_name, vcf[0], vcf[1]] },
            sigprofilermatrixgenerator_interval,
            params.reference,
            params.genome,
            sigprofilermatrixgenerator_reference,
            params.report_s3_dir,
            params.publish_dir,
            params.enable_publish,
            params.disable_publish
        )

        ch_sigprofile_mqc = SIGPROFILERGENERATEMATRIX_WF.out.merged_mutationalcatalog_mqc.collect().ifEmpty([])
        ch_sigprofile_barplot = SIGPROFILERGENERATEMATRIX_WF.out.sbs96_barplot.collect().ifEmpty([])

        if (!params.is_bam) {
            multiqc_finalInput = PRESEQ_WF.out.multiqc_input
                        .combine(ch_sigprofile_mqc)
                        .combine(ch_sigprofile_barplot)
        } else {
            multiqc_finalInput = ch_sigprofile_mqc.combine(ch_sigprofile_barplot)
        }
    } else {
        multiqc_finalInput = PRESEQ_WF.out.multiqc_input
    }

    // DeepVariant subworkflow to call germline variants
    if (params.run_deepvariant) {
        GOOGLE_DEEPVARIANT_WF (
            PRESEQ_WF.out.dedup_bam,
            params.deepvariant_model ? file( params.deepvariant_model, checkIfExists: true ) : [],
            params.population_vcfs   ? file( params.population_vcfs,   checkIfExists: true ) : [],
            params.reference,
            params.intervals         ? file( params.intervals,         checkIfExists: true ) : [],
            params.publish_dir,
            params.enable_publish
        )

        //multiqc_finalInput = multiqc_finalInput.combine(ch_deepvariant_version.collect().ifEmpty([]))
    }

    params_meta = [
            session_id: workflow.sessionId,
            read_length: params.read_length,
            n_reads: params.n_reads,
            genome: params.genome,
            instrument: params.instrument,
            skip_fastqc: params.skip_fastqc,
            skip_kraken: params.skip_kraken,
            skip_qualimap: params.skip_qualimap,
            skip_CNV: params.skip_ginkgo
    ]

    if (!params.skip_ginkgo) {
        params_meta['bin_size'] = params.bin_size
        params_meta['min_ploidy'] = params.min_ploidy
        params_meta['max_ploidy'] = params.max_ploidy
    }

    MULTIQC_WF ( 
                multiqc_finalInput,
                params_meta,
                params.multiqc_config,
                params.publish_dir,
                params.enable_publish
                )
        
    MULTIQC_WF.out.report.ifEmpty{ exit 1, "ERROR: cannot generate any MULTIQC Report." }
}
