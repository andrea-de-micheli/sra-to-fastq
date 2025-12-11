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
    memory '16 GB'
    disk '200 GB'
    
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
