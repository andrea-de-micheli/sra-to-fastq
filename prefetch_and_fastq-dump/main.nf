#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * SRA Prefetch and FASTQ Conversion Pipeline
 * Downloads SRA files and optionally converts to FASTQ format
 */

// Parameters
params.ngc_file = null
params.sra_list_file = null
params.convert_to_fastq = false
params.threads = 4
params.disk_size = 500
params.outdir = "results"

log.info """
    SRA Prefetch and FASTQ Pipeline
    ================================
    ngc_file         : ${params.ngc_file}
    sra_list_file    : ${params.sra_list_file}
    convert_to_fastq : ${params.convert_to_fastq}
    threads          : ${params.threads}
    disk_size        : ${params.disk_size} GB
    outdir           : ${params.outdir}
    """
    .stripIndent()

/*
 * Process: Prefetch SRA files
 */
process PREFETCH_SRA {
    tag "${sra_id}"
    
    container "quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0"
    
    cpus 2
    memory '8 GB'
    disk "${params.disk_size} GB"
    
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    val sra_id
    path ngc_file
    
    output:
    path "${sra_id}", emit: sra
    
    script:
    """
    prefetch ${sra_id} --ngc ${ngc_file} --max-size 0
    """
}

/*
 * Process: Convert SRA to FASTQ
 */
process SRA_TO_FASTQ {
    tag "${sra_file.simpleName}"
    
    container "quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0"
    
    cpus params.threads
    memory '16 GB'
    disk "${params.disk_size} GB"
    
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.fastq.gz"
    
    input:
    path sra_file
    
    output:
    path "*.fastq.gz", emit: fastq
    
    script:
    """
    fasterq-dump \\
        --split-files \\
        --threads ${task.cpus} \\
        --temp . \\
        --outdir . \\
        ${sra_file}
    
    gzip -f *.fastq
    
    # Delete SRA file after successful conversion to save space
    rm -f ${sra_file}
    """
}

/*
 * Main workflow
 */
workflow {
    // Read SRA IDs from file
    sra_list_file = file(params.sra_list_file)
    sra_ch = Channel
        .fromPath(sra_list_file)
        .splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') }
        .ifEmpty { error "No SRA IDs found in ${params.sra_list_file}" }
    
    ngc_file = file(params.ngc_file)
    
    // Run prefetch
    prefetch_out = PREFETCH_SRA(sra_ch, ngc_file)
    
    // Conditionally convert to FASTQ
    if (params.convert_to_fastq) {
        // Find .sra files in the prefetched directories
        // Prefetch creates: SRR_ID/SRR_ID.sra
        sra_files = prefetch_out
            .map { dir -> 
                def sra_id = dir.name
                file("${dir}/${sra_id}.sra")
            }
        
        SRA_TO_FASTQ(sra_files)
    }
}

workflow.onComplete {
    log.info "Pipeline completed: ${workflow.success ? 'OK' : 'Failed'}"
}

