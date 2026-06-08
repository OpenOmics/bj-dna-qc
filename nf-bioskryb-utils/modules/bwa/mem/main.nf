nextflow.enable.dsl=2
params.timestamp = ""

process BWA_MEM {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(reads)
    path fasta_ref
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}_*sorted.bam"), path("${sample_name}_*sorted.bam.bai"), emit: bam
    path("bwa_mem_version.yml"), emit: version
    
    script:
    def memory = "${task.memory.toGiga()-1}G"
    """
    export bwt_max_mem=$memory
 
    bwa mem -M -Y -K 10000000 -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:Illumina" -t $task.cpus '${fasta_ref}/genome.fa' '${reads[0]}' '${reads[1]}' | samtools sort -o '${sample_name}_sorted.bam' -@ $task.cpus
    samtools index '${sample_name}_sorted.bam' -@ $task.cpus
    
    # If switch to BWA-MEM2: \$(echo \$(bwa-mem2 version 2>&1 | sed 's/.* //'))
    cat <<END_VERSIONS > bwa_mem_version.yml
    BWA: \$(echo \$(bwa 2>&1 | sed -n 's/^Version: //p'))
    Samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/ Using.*\$//')
    END_VERSIONS
    """
}

workflow BWA_WF{
    take:
        ch_reads
        ch_reference
        ch_publish_dir
        ch_enable_publish
        
    main:
        BWA_MEM ( 
            ch_reads,
            ch_reference,
            ch_publish_dir,
            ch_enable_publish
        )
                         
    emit:
        bam = BWA_MEM.out.bam
        version = BWA_MEM.out.version
    
}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> def subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], pair] }
    
    BWA_WF ( 
        ch_reads,
        params.reference,
        params.publish_dir,
        params.enable_publish
    )
}