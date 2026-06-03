nextflow.enable.dsl=2
params.timestamp = ""

process SENTIEON_BWA_MEM {
    tag "${sample_name}"
    publishDir "${publish_dir}_${params.timestamp}/${task.process.replaceAll(':', '_')}", enabled:"$enable_publish"

    input:
    tuple val(sample_name), path(reads)
    path fasta_ref
    val(publish_dir)
    val(enable_publish)

    output:
    tuple val(sample_name), path("${sample_name}_*sorted.bam"), path("${sample_name}_*sorted.bam.bai"), emit: bam
    path("sentieon_bwa_mem_version.yml"), emit: version
    
    script:
    def memory = "${task.memory.toGiga()-1}G"
    """
    set +u
    if [ \$LOCAL != "true" ]; then
        . /opt/sentieon/cloud_auth.sh no-op
    else
        export SENTIEON_LICENSE=\$SENTIEON_LICENSE_SERVER
        echo \$SENTIEON_LICENSE
    fi

    export bwt_max_mem=$memory

    # Run Sentieon BWA first and write SAM to disk.
    # This avoids running two Sentieon commands concurrently in a pipe.
    sentieon bwa mem \\
        -M \\
        -Y \\
        -K 10000000 \\
        -R "@RG\\tID:${sample_name}\\tSM:${sample_name}\\tPL:Illumina" \\
        -t $task.cpus \\
        '${fasta_ref}/genome.fa' \\
        '${reads[0]}' \\
        '${reads[1]}' \\
        > '${sample_name}.sam'

    # Then run Sentieon sort after BWA has fully exited.
    sentieon util sort \\
        -r '${fasta_ref}/genome.fa' \\
        -o '${sample_name}_sorted.bam' \\
        -t $task.cpus \\
        --sam2bam \\
        -i '${sample_name}.sam'

    rm -f '${sample_name}.sam' 

    export SENTIEON_VER="202308.01"
    echo Sentieon: \$SENTIEON_VER > sentieon_bwa_mem_version.yml

    """
}

workflow SENTIEON_BWA_WF{
    take:
        ch_reads
        ch_reference
        ch_publish_dir
        ch_enable_publish
        
    main:
        SENTIEON_BWA_MEM ( 
                            ch_reads,
                            ch_reference,
                            ch_publish_dir,
                            ch_enable_publish
                         )
                         
    emit:
        // version = SENTIEON_BWA_MEM.out.version
        bam = SENTIEON_BWA_MEM.out.bam
        version = SENTIEON_BWA_MEM.out.version
    
}

workflow{
    ch_reads = Channel.fromFilePairs( params.reads , size: -1 , checkExists: true )
                            .map { tag, pair -> def subtags = (tag =~ /(.*)_(S\d+)_(L0+\d+)/)[0]; [subtags[1], pair] }
    
    SENTIEON_BWA_WF ( 
                            ch_reads,
                            params.reference,
                            params.publish_dir,
                            params.enable_publish
                         )
}
