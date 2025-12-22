#!/usr/bin/env python3
"""
Preprocess script for SRA Prefetch pipeline.
Extracts uploaded files (NGC key and SRA list) from data URLs and saves them to disk.
"""

import os
import base64
from cirro.helpers.preprocess_dataset import PreprocessDataset

def extract_file_from_data_url(data_url, output_path):
    """Extract file content from data URL and write to disk."""
    if not data_url or not data_url.startswith('data:'):
        raise ValueError(f"Invalid data URL format")
    
    header, encoded = data_url.split(',', 1)
    content = base64.b64decode(encoded)
    
    with open(output_path, 'wb') as f:
        f.write(content)
    
    return os.path.abspath(output_path)

def main():
    ds = PreprocessDataset.from_running()
    
    # Process NGC file
    ngc_data_url = ds.params.get('ngc_file')
    if not ngc_data_url:
        raise ValueError("NGC file is required")
    
    ngc_path = extract_file_from_data_url(ngc_data_url, "ngc_key.ngc")
    ds.remove_param("ngc_file")
    ds.add_param("ngc_file", ngc_path, overwrite=True)
    
    # Process SRA list file
    sra_list_data_url = ds.params.get('sra_list_file')
    if not sra_list_data_url:
        raise ValueError("SRA list file is required")
    
    sra_list_path = extract_file_from_data_url(sra_list_data_url, "sra_list.txt")
    ds.remove_param("sra_list_file")
    ds.add_param("sra_list_file", sra_list_path, overwrite=True)
    
    # Count SRA IDs for logging
    with open(sra_list_path, 'r') as f:
        sra_count = len([line.strip() for line in f if line.strip() and not line.strip().startswith('#')])
    
    print(f"Preprocessing complete: {sra_count} SRA ID(s) found")

if __name__ == "__main__":
    main()

