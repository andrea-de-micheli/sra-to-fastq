#!/usr/bin/env python3
"""
Preprocess script for SRA to FASTQ pipeline.

This script runs BEFORE the Nextflow pipeline starts.
It can be used to:
- Validate inputs
- Create samplesheets
- Modify parameters
"""
import base64
import os
from cirro.helpers.preprocess_dataset import PreprocessDataset

def main():
    # Initialize the dataset helper
    ds = PreprocessDataset.from_running()
    
    # Log what we're working with
    ds.logger.info(f"Parameters: {ds.params}")
    
    # Handle sample sheet upload if provided
    samplesheet_file = ds.params.get('samplesheet_file')
    if samplesheet_file and samplesheet_file.startswith('data:'):
        try:
            # Decode data URL (format: data:text/csv;base64,<content>)
            _, encoded = samplesheet_file.split(',', 1)
            content = base64.b64decode(encoded)
            
            # Save as samplesheet.csv in dataset root
            samplesheet_path = os.path.join(ds.dataset_root, "samplesheet.csv")
            with open(samplesheet_path, 'wb') as f:
                f.write(content)
            ds.logger.info(f"Sample sheet saved to: {samplesheet_path}")
        except Exception as e:
            ds.logger.error(f"Error saving sample sheet: {e}")
    
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

if __name__ == "__main__":
    main()