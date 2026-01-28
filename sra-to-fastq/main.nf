#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * SRA to FASTQ Conversion Pipeline
 * Converts .sra files to compressed FASTQ format using fasterq-dump
 */

// Parameters (populated by Cirro via process-input.json)
params.input_dir = null      // S3 path to input dataset
params.threads = 4           // CPU threads for fasterq-dump
params.outdir = "results"    // Output directory
params.skip_sra_ids = null   // Comma-separated list of SRA IDs to skip

// Log parameters
log.info """
    SRA to FASTQ Pipeline
    =====================
    input_dir : ${params.input_dir}
    threads   : ${params.threads}
    outdir    : ${params.outdir}
    """
    .stripIndent()

/*
 * Process: Extract FASTQ from SRA
 */
process EXTRACT_FASTQ {
    tag "${sra_file.simpleName}"
    
    // Use BioContainers image for sra-tools
    container "quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0"
    
    // Resource allocation
    cpus params.threads
    memory '32 GB'
    
    input:
    path sra_file
    
    output:
    // One or more FASTQ files per SRA (e.g. *_1.fastq, *_2.fastq)
    path "*.fastq", emit: fastq_files
    
    script:
    """
    fasterq-dump \\
        --split-files \\
        --threads ${task.cpus} \\
        --disk-limit-tmp 0 \\
        --disk-limit 0 \\
        --temp . \\
        --outdir . \\
        ${sra_file}
    """
}

/*
 * Process: Compress FASTQ with pigz
 */
process COMPRESS_FASTQ {
    tag "${fastq.simpleName}"
    
    // Use a container that has pigz
    container "quay.io/biocontainers/pigz:2.8--h2797004_0"
    
    cpus params.threads
    memory '8 GB'
    
    // Publish compressed outputs to results directory
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.fastq.gz"
    
    input:
    path fastq
    
    output:
    path "*.fastq.gz"
    
    script:
    """
    pigz -p ${task.cpus} ${fastq}
    """
}

/*
 * Main workflow
 */
workflow {
    // Find all .sra files recursively in the input directory
    sra_ch = Channel
        .fromPath("${params.input_dir}/**/*.sra", checkIfExists: true)
        .ifEmpty { error "No .sra files found in ${params.input_dir}" }
    
    // Filter out skipped SRA files if skip_sra_ids is provided
    // Supports both comma-separated and newline-separated lists
    def skip_list = []
    if (params.skip_sra_ids) {
        // Split by both comma and newline, then filter out empty strings
        skip_list = params.skip_sra_ids
            .split(/[,\n\r]+/)
            .collect { it.trim() }
            .findAll { it.length() > 0 }
    }
    
    if (skip_list) {
        log.info "Skipping ${skip_list.size()} SRA IDs: ${skip_list.join(', ')}"
        sra_ch = sra_ch.filter { sra_file ->
            def sra_id = sra_file.getName().replace('.sra', '')
            def should_skip = skip_list.contains(sra_id)
            if (should_skip) {
                log.info "Skipping SRA file: ${sra_file}"
            }
            !should_skip
        }
    }
    
    // Log what we found
    sra_ch.view { "Found SRA file: ${it}" }
    
    /*
     * Stage 1: Extract FASTQ
     */
    EXTRACT_FASTQ(sra_ch)
    def fastq_ch = EXTRACT_FASTQ.out.fastq_files
    
    // Flatten to handle multiple FASTQ files per SRA
    fastq_ch
        .flatten()
        .set { individual_fastq }
    
    /*
     * Stage 2: Compress each FASTQ in parallel
     */
    COMPRESS_FASTQ(individual_fastq)
}

// On completion
workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'Failed'}"
}
