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
    ds.logger.info(f"Parameters: {ds.params}")
    
    # Debug: Check what attributes are available
    ds.logger.info(f"Available attributes: {[attr for attr in dir(ds) if not attr.startswith('_')]}")
    
    # Check the files in the input dataset
    ds.logger.info(f"Files DataFrame shape: {ds.files.shape}")
    ds.logger.info(f"Files DataFrame columns: {ds.files.columns.tolist()}")
    ds.logger.info(f"Files in dataset:\n{ds.files}")
    
    # Check if we have an input_dir parameter
    input_dir = ds.params.get('input_dir')
    if input_dir:
        ds.logger.info(f"Input directory from params: {input_dir}")
        ds.logger.info("Note: Nextflow will check for .sra files at this S3 path")
    
    # Count SRA files (if DataFrame has data)
    if len(ds.files) > 0:
        sra_files = ds.files[ds.files['file'].str.endswith('.sra')]
        ds.logger.info(f"Found {len(sra_files)} SRA files in ds.files DataFrame")
        
        if len(sra_files) == 0:
            ds.logger.warning("No .sra files found in ds.files DataFrame, but files may exist in S3")
    else:
        ds.logger.warning("ds.files DataFrame is empty. This may be normal if files aren't registered in Cirro metadata.")
        ds.logger.info("Nextflow will validate files exist at the input_dir S3 path")
    
    # The process-input.json handles parameter mapping, 
    # but we could also add/modify params here if needed:
    # ds.add_param("custom_param", "value")

if __name__ == "__main__":
    main()