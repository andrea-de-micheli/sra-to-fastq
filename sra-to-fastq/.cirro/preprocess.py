#!/usr/bin/env python3
"""
Preprocess script for SRA to FASTQ pipeline.

This script runs BEFORE the Nextflow pipeline starts.
It can be used to:
- Validate inputs
- Create samplesheets
- Modify parameters

For this simple pipeline, we just validate that SRA files exist.
"""
from cirro.helpers.preprocess_dataset import PreprocessDataset

def main():
    # Initialize the dataset helper
    ds = PreprocessDataset.from_running()
    
    # Log what we're working with
    ds.logger.info(f"Dataset ID: {ds.dataset['id']}")
    ds.logger.info(f"Dataset S3 path: {ds.dataset['s3']}")
    ds.logger.info(f"Parameters: {ds.params}")
    
    # Check the files in the input dataset
    ds.logger.info(f"Files in dataset:\n{ds.files}")
    
    # Count SRA files
    sra_files = ds.files[ds.files['name'].str.endswith('.sra')]
    ds.logger.info(f"Found {len(sra_files)} SRA files to convert")
    
    if len(sra_files) == 0:
        raise ValueError("No .sra files found in input dataset!")
    
    # The process-input.json handles parameter mapping, 
    # but we could also add/modify params here if needed:
    # ds.add_param("custom_param", "value")

if __name__ == "__main__":
    main()
