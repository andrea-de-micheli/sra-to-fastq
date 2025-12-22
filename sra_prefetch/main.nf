#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * SRA Prefetch Pipeline
 * Downloads SRA files using prefetch with NGC key authentication
 */

// Parameters (populated by Cirro via process-input.json)
params.ngc_file = null      // Path to NGC key file
params.sra_list_file = null // Path to file containing list of SRA IDs (one per line)
params.max_file_size = 500  // Max expected file size in GB (default: 500 GB)
params.outdir = "results"   // Output directory

// Log parameters
log.info """
    SRA Prefetch Pipeline
    =====================
    ngc_file      : ${params.ngc_file}
    sra_list_file : ${params.sra_list_file}
    max_file_size : ${params.max_file_size} GB
    outdir        : ${params.outdir}
    """
    .stripIndent()

/*
 * Process: Prefetch SRA files
 */
process PREFETCH_SRA {
    tag "${sra_id}"
    
    // Use BioContainers image for sra-tools
    container "quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0"
    
    // Resource allocation
    cpus 2
    memory '8 GB'
    disk "${params.max_file_size} GB"
    
    // Publish outputs to results directory
    // Prefetch creates a directory for each SRA, so we copy the entire directory
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    val sra_id
    path ngc_file
    
    output:
    path "${sra_id}", emit: sra
    
    script:
    """
    # Run prefetch with NGC key (--max-size 0 removes size limit)
    prefetch ${sra_id} --ngc ${ngc_file} --max-size 0
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
    
    // Load NGC file
    ngc_file = file(params.ngc_file)
    
    // Log what we found
    sra_ch.view { "Prefetching SRA: ${it}" }
    
    // Run prefetch for each SRA
    PREFETCH_SRA(sra_ch, ngc_file)
}

// On completion
workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'Failed'}"
}

