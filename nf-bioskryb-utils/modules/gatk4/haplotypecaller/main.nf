nextflow.enable.dsl=2
params.timestamp = ""

process GATK4_HAPLOTYPECALLER {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/secondary_analyses/variant_calls_haplotypecaller", enabled:"$enable_publish", pattern: "*_haplotypecaller.vcf.gz*"

    input:
    tuple val(sample_name), path(bam), path(bai)
    path(fasta_ref)
    path(interval)
    path(dbsnp)
    path(dbsnp_index)
    val(ploidy)
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("*_haplotypecaller.vcf.gz*"), emit: vcf
    path("gatk4_haplotypecaller_version.yml"), emit: version

    script:
    def dbsnp_param = dbsnp ? "--dbsnp ${dbsnp}" : ""
    def interval_param = interval ? "--intervals ${interval}" : ""
    def avail_mem = (task.memory.mega * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}M" HaplotypeCaller \\
        --reference ${fasta_ref}/genome.fa \\
        --input ${bam} \\
        --output ${sample_name}_haplotypecaller.vcf.gz \\
        --native-pair-hmm-threads ${task.cpus} \\
        --sample-ploidy ${ploidy} \\
        ${dbsnp_param} \\
        ${interval_param}
    
    export GATK4_VER=\$(echo \$(gatk --version 2>&1) | sed -e 's/.*(GATK) //; s/Version: //g')
    echo GATK4: \$GATK4_VER > gatk4_haplotypecaller_version.yml
    """
}

