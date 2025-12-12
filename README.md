# Cirro Custom Nextflow Pipelines Guide

A guide for creating and deploying Nextflow pipelines on Cirro.

For detailed Cirro documentation, see the [Cirro Documentation](https://docs.cirro.bio/).

## How Cirro Pipelines Work

When you run a pipeline in Cirro, here's what happens:

```
1. User fills out web form (process-form.json)
           ↓
2. Form values mapped to parameters (process-input.json)
           ↓
3. Preprocess script runs (preprocess.py) [OPTIONAL]
   - Can modify parameters
   - Can create samplesheets
   - Has access to input dataset info
           ↓
4. Nextflow workflow executes (main.nf)
           ↓
5. Output files saved to new dataset
```

## Repository Structure

Each Nextflow pipelines needs to be in its own dedicated directory and needs a `.cirro/` folder with configuration files:

```
your-nf-pipeline/
├── main.nf                      # Your Nextflow pipeline
├── nextflow.config              # Nextflow configuration
└── .cirro/
    ├── process-form.json        # REQUIRED: Defines the web form UI
    ├── process-input.json       # REQUIRED: Maps form → Nextflow params
    ├── process-output.json      # OPTIONAL: For dashboard visualizations
    ├── process-compute.config   # OPTIONAL: Override Nextflow settings
    └── preprocess.py            # OPTIONAL: Pre-workflow Python logic
```

## Configuration Files

### 1. process-form.json (REQUIRED)

Defines what the user sees in the web interface. Uses JSON Schema format. This file controls the UI form that users interact with when running your pipeline.

```json
{
    "ui": {},
    "form": {
        "title": "My Pipeline",
        "type": "object",
        "properties": {
            "my_param": {
                "title": "Human Readable Name",
                "description": "Help text for user",
                "type": "integer",
                "default": 4
            }
        }
    }
}
```

For more information on form configuration, see the [Cirro Custom Pipelines documentation](https://docs.cirro.bio/).

### 2. process-input.json (REQUIRED)

Maps form values to Nextflow parameters. This is how user inputs get passed to your pipeline.

```json
{
    "my_param": "$.dataset.params.my_param",
    "input_dir": "$.inputs[*].dataPath"
}
```

Key JSONPath expressions:
- `$.dataset.params.X` → User form values
- `$.inputs[*].dataPath` → S3 path(s) to input dataset(s)
- `$.dataset.dataPath` → S3 path for output dataset

### 3. preprocess.py (OPTIONAL)

Runs BEFORE Nextflow execution. Use it to:
- Create samplesheets dynamically
- Modify parameters based on input data
- Validate inputs
- Generate configuration files

```python
from cirro.helpers.preprocess_dataset import PreprocessDataset

ds = PreprocessDataset.from_running()

# Access dataset info:
ds.files        # DataFrame of files in input dataset
ds.samplesheet  # Sample metadata from samplesheet.csv (if present)
ds.params       # User form parameters
ds.dataset      # Full dataset object (includes S3 paths)

# Modify parameters:
ds.add_param("input_files", "/path/to/files")
ds.remove_param("unwanted_param")

# Write files (e.g., samplesheet for Nextflow):
df.to_csv("samplesheet.csv", index=False)
```

### 4. process-compute.config (OPTIONAL)

Override Nextflow settings for compute resources (memory, containers, etc.). This file is merged with your base `nextflow.config`.

```groovy
params {
    my_container = "quay.io/biocontainers/tool:version"
}

process {
    withName: 'MY_PROCESS' {
        memory = '32 GB'
        cpus = 8
    }
}
```

### 5. process-output.json (OPTIONAL)

Defines dashboard visualizations for pipeline outputs. See the [Cirro Custom Pipelines documentation](https://docs.cirro.bio/) for details on output visualization configuration.

## Nextflow Pipeline Basics

Your `main.nf` file should follow standard Nextflow conventions:

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters come from process-input.json
params.input_dir = null   // S3 path to input dataset
params.threads = 4        // From user form

// A "process" is a single computational step
process MY_STEP {
    container 'quay.io/biocontainers/tool:version'
    cpus params.threads
    memory '16 GB'
    
    publishDir "results", mode: 'copy'  // Where outputs go
    
    input:
    path input_file
    
    output:
    path "*.output"
    
    script:
    """
    my_tool ${input_file}
    """
}

// The workflow connects processes
workflow {
    // Create a channel (stream of data)
    input_ch = Channel.fromPath("${params.input_dir}/**/*.input")
    
    // Run the process on each item
    MY_STEP(input_ch)
}
```

### Key Nextflow Concepts

1. **Channels**: Streams of data that flow through your pipeline
2. **Processes**: Individual computational steps (each runs in its own container)
3. **publishDir**: Where output files are saved (relative to the output dataset)
4. **Containers**: Each process runs in a Docker container

For comprehensive Nextflow documentation, see [Nextflow Documentation](https://www.nextflow.io/docs/latest/).

## Finding Containers

Use BioContainers for pre-built bioinformatics tools:
- Browse: https://biocontainers.pro/
- Direct: `quay.io/biocontainers/<tool>:<version>`