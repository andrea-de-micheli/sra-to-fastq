# H5AD Single-Cell QC (Cirro + Nextflow)

This pipeline scans an input dataset for all `.h5ad` files, runs single-cell QC preprocessing on each file, and writes results into one subdirectory per input file.

## Inputs

- Input dataset path (mapped to `input_dir`)
- Cirro form parameters:
  - `min_genes` (default: `200`)
  - `min_counts` (default: `500`)
  - `n_mads` (default: `3.0`)
  - `mt_cutoff` (default: `null` / MAD-based MT threshold)
  - `expected_doublet_rate` (default: `0.06`)
  - `skip_doublets` (default: `false`)

## Outputs

For each input file `<sample>.h5ad`, the pipeline writes:

- `<outdir>/<sample>/<sample>.cleaned.h5ad`
- `<outdir>/<sample>/<sample>.qc_stats.txt`
- `<outdir>/<sample>/<sample>.qc_plots.png`

## Implementation Notes

- Workflow entrypoint: `main.nf`
- Helper script: `bin/single_cell_preprocessing_pipeline.py`
- Container: `getwilds/scanpy:1.10.2`
- Runtime dependency install in process: `scrublet==0.2.3`
