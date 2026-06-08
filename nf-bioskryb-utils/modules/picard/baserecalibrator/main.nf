nextflow.enable.dsl=2
params.timestamp = ""

process GATK_BQSR {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/alignment", enabled:"$enable_publish"

    input:
    val(genome)
    tuple val(sample_name), path(bam), path(bai)
    path fasta_ref
    path dbsnp
    path dbsnp_index
    path mills
    path mills_index
    path onekg
    path onekg_index
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}_recal_data.table"), emit: recal_table
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai"), path("${sample_name}_recal_data.table"), emit: bam_recal_table
    path("gatk4_bsqr_version.yml"), emit: version

    script:
    """
    if [[ ${genome} =~ .*GRCh3* ]];
    then
        # Generate first table
        gatk BaseRecalibrator \
            -R ${fasta_ref}/genome.fa \
            -I ${bam} \
            --known-sites ${dbsnp} \
            --known-sites ${mills} \
            --known-sites ${onekg} \
            -O ${sample_name}_recal_data.table

        # Apply recalibration
        gatk ApplyBQSR \
            -R ${fasta_ref}/genome.fa \
            -I ${bam} \
            --bqsr-recal-file ${sample_name}_recal_data.table \
            -O ${sample_name}_recal.bam
        
        # Generate post-recalibration table
         gatk BaseRecalibrator \
            -R ${fasta_ref}/genome.fa \
            -I ${sample_name}_recal.bam \
            --known-sites ${dbsnp} \
            --known-sites ${mills} \
            --known-sites ${onekg} \
            -O ${sample_name}_recal_data.table.post

    else
        # Generate first table
        gatk BaseRecalibrator \
            -R ${fasta_ref}/genome.fa \
            -I ${bam} \
            -O ${sample_name}_recal_data.table

        # Apply recalibration
        gatk ApplyBQSR \
            -R ${fasta_ref}/genome.fa \
            -I ${bam} \
            --bqsr-recal-file ${sample_name}_recal_data.table \
            -O ${sample_name}_recal.bam
        
        # Generate post-recalibration table
         gatk BaseRecalibrator \
            -R ${fasta_ref}/genome.fa \
            -I ${sample_name}_recal.bam \
            -O ${sample_name}_recal_data.table.post
            
    fi
        
    # Plot BQSR outcomes
    gatk AnalyzeCovariates \
        -before ${sample_name}_recal_data.table \
        -after ${sample_name}_recal_data.table.post \
        -plots ${sample_name}_recal_plots.pdf
    
    export GATK4_VER=\$(echo \$(gatk --version 2>&1) | sed -e 's/.*(GATK) //; s/Version: //g')
    echo GATK4: \$GATK4_VER > gatk4_bsqr_version.yml
    """
}