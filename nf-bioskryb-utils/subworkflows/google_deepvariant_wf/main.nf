/*
========================================================================================
 Google DeepVariant — shareable Nextflow subworkflow
----------------------------------------------------------------------------------------
 Three-stage DeepVariant implementation:
   1. DEEPVARIANT_MAKE_EXAMPLES_ONLY  — CPU, shards examples across task.cpus
   2. DEEPVARIANT_CALL_VARIANTS       — typically GPU, runs the model
   3. DEEPVARIANT_POSTPROCESS         — CPU, writes the final VCF

 Self-contained: no includeConfig to internal shared configs, no private
 container registries, no private reference bucket references.
 Uses the public `google/deepvariant` images from Docker Hub.
========================================================================================
*/

nextflow.enable.dsl = 2
params.timestamp = ""


process DEEPVARIANT_MAKE_EXAMPLES_ONLY {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(bam), path(bai)
    path(deepvariant_model)
    path(population_vcfs)
    path(reference)
    path(regions)

    output:
    tuple val("${sample_name}"), path("${sample_name}.tfrecord-*-of-*.gz"), val(task.cpus), emit: example_tfrecords

    script:
    def regions_arg        = regions        ? "--regions ${regions}" : ""
    def checkpoint_arg     = deepvariant_model ? "--checkpoint ${deepvariant_model}/${params.checkpoint_filename}" : ""
    def pop_vcfs_setup     = population_vcfs ? """population_vcfs_list=\$(ls ${population_vcfs}/*.vcf.gz 2>/dev/null | tr '\\n' ',' | sed 's/,\$//' || echo "")""" : ""
    def pop_vcfs_arg       = population_vcfs ? '--population_vcfs="${population_vcfs_list}"' : ""
    """
    set -e
    export DV_BIN_PATH=/opt/deepvariant/bin
    mkdir -p scratch_${sample_name}
    export TMPDIR=\$PWD/scratch_${sample_name}

    ${pop_vcfs_setup}

    seq 0 ${task.cpus - 1} | \\
    parallel -q --halt 2 --line-buffer \\
    /opt/deepvariant/bin/make_examples \\
      --mode calling \\
      --ref ${reference}/${params.reference_fasta_name} \\
      --reads ${bam} \\
      --examples "${sample_name}.tfrecord@${task.cpus}.gz" \\
      ${checkpoint_arg} \\
      ${pop_vcfs_arg} \\
      ${regions_arg} \\
      --task {}
    """
}


process DEEPVARIANT_CALL_VARIANTS {
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(example_tfrecords), val(shards)
    path(deepvariant_model)

    output:
    tuple val(sample_name), path("${sample_name}_variants_output*tfrecord.gz"), path("variants_shards.txt"), emit: variants_output

    script:
    def checkpoint_arg = deepvariant_model ? "--checkpoint ${deepvariant_model}/${params.checkpoint_filename}" : ""
    """
    set -e

    /opt/deepvariant/bin/call_variants \\
      --examples=${sample_name}.tfrecord@${shards}.gz \\
      --outfile=${sample_name}_variants_output.tfrecord.gz \\
      ${checkpoint_arg} \\
      --batch_size=${params.call_variants_batch_size}

    # Count files produced (for postprocess_variants' shard pattern)
    ls -1 ${sample_name}_variants_output*.tfrecord.gz | wc -l > variants_shards.txt
    """
}


process DEEPVARIANT_POSTPROCESS {
    tag "${sample_name}"
    publishDir  "${params.publish_dir}_${params.timestamp}/deepvariant/", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(variants_output), path(variants_shards_file)
    path(reference)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}_deepvariant.vcf.gz*"), emit: vcf
    path("deepvariant_version.yml"),                                    emit: version

    script:
    """
    set -e

    # Read the actual shard count written by DEEPVARIANT_CALL_VARIANTS
    variants_output_shards=\$(cat ${variants_shards_file} | tr -d '[:space:]')

    if [ -z "\$variants_output_shards" ] || [ "\$variants_output_shards" -lt 1 ]; then
        echo "Error: Invalid shard count: \$variants_output_shards" >&2
        exit 1
    fi

    /opt/deepvariant/bin/postprocess_variants \\
      --ref=${reference}/${params.reference_fasta_name} \\
      --infile=${sample_name}_variants_output@\${variants_output_shards}.tfrecord.gz \\
      --outfile="${sample_name}_deepvariant.vcf.gz" \\
      --cpus=${task.cpus}

    # Record the DeepVariant image/version for provenance
    echo "DeepVariant: ${params.deepvariant_version}" > deepvariant_version.yml
    """
}


workflow GOOGLE_DEEPVARIANT_WF {

    take:
        ch_bam              // channel: [ sample_name, bam, bai ]
        ch_deepvariant_model   // value:   directory with model files (or [] for container default)
        ch_population_vcfs     // value:   directory with per-chromosome VCFs (or [] to skip)
        ch_reference           // value:   directory with reference FASTA + .fai
        ch_regions             // value:   BED file (or [] for whole-genome)
        ch_publish_dir
        ch_enable_publish

    main:
        DEEPVARIANT_MAKE_EXAMPLES_ONLY(
            ch_bam,
            ch_deepvariant_model,
            ch_population_vcfs,
            ch_reference,
            ch_regions
        )

        DEEPVARIANT_CALL_VARIANTS(
            DEEPVARIANT_MAKE_EXAMPLES_ONLY.out.example_tfrecords,
            ch_deepvariant_model
        )

        DEEPVARIANT_POSTPROCESS(
            DEEPVARIANT_CALL_VARIANTS.out.variants_output,
            ch_reference,
            ch_publish_dir,
            ch_enable_publish
        )

    emit:
        vcf     = DEEPVARIANT_POSTPROCESS.out.vcf
        version = DEEPVARIANT_POSTPROCESS.out.version
}


workflow {

    // ------------------------------------------------------------------
    // Parameter validation
    // ------------------------------------------------------------------
    if ( !params.input_csv ) { error "Missing required parameter: --input_csv" }
    if ( !params.reference ) { error "Missing required parameter: --reference (directory containing the reference fasta)" }

    // ------------------------------------------------------------------
    // Build input channels
    //   Expected CSV columns: sample_name,bam[,bai]
    //   (bai auto-derived from <bam>.bai when not provided)
    // ------------------------------------------------------------------
    ch_bam = Channel.fromPath( params.input_csv, checkIfExists: true )
        .splitCsv( header: true, sep: ',' )
        .map { row ->
            def bam = file( row.bam, checkIfExists: true )
            def bai = row.bai ? file( row.bai, checkIfExists: true ) : file( "${row.bam}.bai", checkIfExists: true )
            tuple( row.sample_name, bam, bai )
        }

    ch_reference         = file( params.reference, checkIfExists: true )
    ch_deepvariant_model = params.deepvariant_model ? file( params.deepvariant_model, checkIfExists: true ) : []
    ch_population_vcfs   = params.population_vcfs   ? file( params.population_vcfs,   checkIfExists: true ) : []
    ch_regions           = params.regions           ? file( params.regions,           checkIfExists: true ) : []

    GOOGLE_DEEPVARIANT_WF(
        ch_bam,
        ch_deepvariant_model,
        ch_population_vcfs,
        ch_reference,
        ch_regions
    )

    GOOGLE_DEEPVARIANT_WF.out.vcf.view()
}
