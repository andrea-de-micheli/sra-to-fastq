#!/usr/bin/env python3
"""
Preprocess script for SRA Prefetch pipeline.
Handles NGC file upload and SRA list processing.
"""

import os
import base64
from pathlib import Path
from cirro.helpers.preprocess_dataset import PreprocessDataset

def main():
    ds = PreprocessDataset.from_running()
    
    # Get NGC file from params (it comes as a data URL)
    ngc_data_url = ds.params.get('ngc_file')
    if not ngc_data_url:
        raise ValueError("NGC file is required")
    
    # Extract the file content from data URL
    # Format: data:application/octet-stream;base64,<base64_content>
    if ngc_data_url.startswith('data:'):
        # Parse data URL
        header, encoded = ngc_data_url.split(',', 1)
        # Extract content type if needed
        content_type = header.split(';')[0].split(':')[1] if ':' in header else None
        
        # Decode base64 content
        ngc_content = base64.b64decode(encoded)
        
        # Write NGC file
        ngc_file_path = "ngc_key.ngc"
        with open(ngc_file_path, 'wb') as f:
            f.write(ngc_content)
        
        # Add as parameter (absolute path)
        ngc_abs_path = os.path.abspath(ngc_file_path)
        ds.add_param("ngc_file", ngc_abs_path)
        ngc_file_display = ngc_file_path
    else:
        # Assume it's already a file path
        ngc_abs_path = os.path.abspath(ngc_data_url) if not os.path.isabs(ngc_data_url) else ngc_data_url
        ds.add_param("ngc_file", ngc_abs_path)
        ngc_file_display = ngc_data_url
    
    # Get SRA list from params
    sra_list = ds.params.get('sra_list', '')
    if not sra_list or not sra_list.strip():
        raise ValueError("SRA list is required")
    
    # Write SRA list to file (one per line)
    sra_list_path = "sra_list.txt"
    with open(sra_list_path, 'w') as f:
        # Split by newlines and write each non-empty line
        for line in sra_list.strip().split('\n'):
            line = line.strip()
            if line and not line.startswith('#'):
                f.write(line + '\n')
    
    # Add as parameter (absolute path)
    sra_list_abs_path = os.path.abspath(sra_list_path)
    ds.add_param("sra_list_file", sra_list_abs_path)
    
    # Log summary
    sra_count = len([l for l in sra_list.strip().split('\n') if l.strip() and not l.strip().startswith('#')])
    print(f"Preprocessing complete:")
    print(f"  - NGC file: {ngc_file_display}")
    print(f"  - SRA IDs to prefetch: {sra_count}")

if __name__ == "__main__":
    main()

