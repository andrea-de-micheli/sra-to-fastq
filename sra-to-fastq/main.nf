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
 * Process: Convert SRA to FASTQ
 */
process SRA_TO_FASTQ {
    tag "${sra_file.simpleName}"
    
    // Use BioContainers image for sra-tools
    container "quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0"
    
    // Resource allocation
    cpus params.threads
    memory '32 GB'
    disk '2500 GB'
    
    // Publish outputs to results directory
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.fastq.gz"
    
    input:
    path sra_file
    
    output:
    path "*.fastq.gz", emit: fastq
    
    script:
    def prefix = sra_file.simpleName
    """
    # Run fasterq-dump with split-files for paired-end support
    fasterq-dump \\
        --split-files \\
        --threads ${task.cpus} \\
        --temp . \\
        --outdir . \\
        ${sra_file}
    
    # Compress all resulting FASTQ files
    gzip -f *.fastq
    
    # Rename files to use clean sample names if needed
    # (fasterq-dump already uses the SRA ID as prefix)
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
    
    // Run the conversion
    SRA_TO_FASTQ(sra_ch)
}

// On completion
workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'Failed'}"
}
