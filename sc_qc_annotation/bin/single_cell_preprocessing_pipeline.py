#!/usr/bin/env python3
"""
Single-cell RNA-seq preprocessing script.

Loads a raw .h5ad file, computes QC metrics (total UMI counts, genes detected,
mitochondrial %, cell cycle scores), generates QC distribution plots with
threshold lines, filters low-quality cells by percentile-based cutoffs, removes
doublets with Scrublet, and saves the cleaned AnnData object.

Dependencies (pip install):
    scanpy anndata scrublet matplotlib numpy

Usage:
    python preprocess_sc.py --input raw.h5ad --output cleaned.h5ad

Run with --help for the full list of options.
"""

import argparse

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scrublet as scr

# ---------------------------------------------------------------------------
# Tirosh 2015 cell-cycle gene sets (same lists that Seurat's cc.genes uses).
# Two separate lists for S-phase and G2/M-phase genes.
# ---------------------------------------------------------------------------

S_GENES_HUMAN = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

G2M_GENES_HUMAN = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
    "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
    "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_adata(path: str) -> ad.AnnData:
    """Load an h5ad file and ensure .X holds raw integer counts."""
    print(f"Loading {path} ...")
    adata = sc.read_h5ad(path)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # If the object has a .raw slot and .X looks normalized, swap in raw counts.
    if adata.raw is not None:
        x_max = np.max(adata.X[:100].toarray() if hasattr(adata.X, "toarray") else adata.X[:100])
        if x_max < 50:
            print("  .X looks normalized (max in first 100 cells < 50). Swapping in .raw.X as counts.")
            adata = adata.raw.to_adata()

    # Remove cells with zero total counts (empty droplets / barcodes with no RNA).
    total_counts = np.array(adata.X.sum(axis=1)).flatten()
    n_zero = (total_counts == 0).sum()
    if n_zero > 0:
        print(f"  Removing {n_zero} cells with zero total counts.")
        adata = adata[total_counts > 0].copy()

    return adata


def compute_qc_metrics(adata: ad.AnnData) -> None:
    """Compute per-cell QC metrics in-place: total counts, gene counts, MT%."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_mt = adata.var["mt"].sum()
    print(f"  Found {n_mt} mitochondrial genes (prefix 'MT-')")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    # This populates:
    #   adata.obs["total_counts"]         ~ nCount_RNA
    #   adata.obs["n_genes_by_counts"]    ~ nFeature_RNA
    #   adata.obs["pct_counts_mt"]        ~ percent.mt


def score_cell_cycle(adata: ad.AnnData) -> None:
    """Score each cell for S and G2M phase on a temporary normalized copy,
    then write scores back into the original adata.obs."""
    s_genes = S_GENES_HUMAN
    g2m_genes = G2M_GENES_HUMAN

    # Work on a copy so we don't alter the raw counts in .X
    adata_tmp = adata.copy()
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)

    # Keep only genes present in the dataset
    s_found = [g for g in s_genes if g in adata_tmp.var_names]
    g2m_found = [g for g in g2m_genes if g in adata_tmp.var_names]
    print(f"  Cell-cycle genes found: {len(s_found)}/{len(s_genes)} S-phase, "
          f"{len(g2m_found)}/{len(g2m_genes)} G2M-phase")

    if len(s_found) < 5 or len(g2m_found) < 5:
        print("  WARNING: Too few cell-cycle genes found. Skipping scoring.")
        adata.obs["S_score"] = np.nan
        adata.obs["G2M_score"] = np.nan
        adata.obs["phase"] = "Unknown"
        adata.obs["CC_difference"] = np.nan
        return

    sc.tl.score_genes_cell_cycle(adata_tmp, s_genes=s_found, g2m_genes=g2m_found)

    adata.obs["S_score"] = adata_tmp.obs["S_score"].values
    adata.obs["G2M_score"] = adata_tmp.obs["G2M_score"].values
    adata.obs["phase"] = adata_tmp.obs["phase"].values
    adata.obs["CC_difference"] = adata.obs["S_score"] - adata.obs["G2M_score"]


def mad_thresholds(values: np.ndarray, n_mads: float = 3) -> tuple[float, float]:
    """Compute median ± n_mads * MAD on log1p-transformed values.

    Returns (lo, hi) thresholds in the original (non-log) scale.
    """
    log_vals = np.log1p(values)
    median = np.median(log_vals)
    mad = np.median(np.abs(log_vals - median))
    lo = np.expm1(max(0, median - n_mads * mad))
    hi = np.expm1(median + n_mads * mad)
    return lo, hi


def compute_thresholds(adata: ad.AnnData, args) -> dict:
    """Compute QC cutoffs using MAD-based outlier detection."""
    counts = adata.obs["total_counts"].values
    genes = adata.obs["n_genes_by_counts"].values
    mt = adata.obs["pct_counts_mt"].values
    n_mads = args.n_mads

    # MAD-based thresholds for UMI and gene counts (log-space)
    counts_lo, counts_hi = mad_thresholds(counts, n_mads)
    genes_lo, genes_hi = mad_thresholds(genes, n_mads)

    # Apply hard floors
    counts_lo = max(args.min_counts, counts_lo)
    genes_lo = max(args.min_genes, genes_lo)

    # MT%: use absolute cutoff if provided, otherwise MAD-based (non-log)
    if args.mt_cutoff is not None:
        mt_hi = args.mt_cutoff
        mt_method = "absolute cutoff"
    else:
        median_mt = np.median(mt)
        mad_mt = np.median(np.abs(mt - median_mt))
        mt_hi = median_mt + n_mads * mad_mt
        mt_method = f"MAD-based, n_mads={n_mads}"

    thresholds = {
        "counts_lo": counts_lo,
        "counts_hi": counts_hi,
        "genes_lo": genes_lo,
        "genes_hi": genes_hi,
        "mt_hi": mt_hi,
    }

    print(f"\n--- QC Thresholds (MAD-based, n_mads={n_mads}) ---")
    print(f"  UMI counts : {thresholds['counts_lo']:.0f} – {thresholds['counts_hi']:.0f}  "
          f"(floor={args.min_counts})")
    print(f"  Gene counts: {thresholds['genes_lo']:.0f} – {thresholds['genes_hi']:.0f}  "
          f"(floor={args.min_genes})")
    print(f"  MT %        : ≤ {thresholds['mt_hi']:.2f}%  "
          f"({mt_method})")

    return thresholds


def plot_qc(adata: ad.AnnData, thresholds: dict, doublet_scores: np.ndarray | None,
            doublet_threshold: float | None, output_path: str) -> None:
    """Generate a multi-panel QC figure and save to disk."""
    n_panels = 6 if doublet_scores is not None else 5
    n_cols = 3
    n_rows = 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 10))
    axes = axes.flatten()

    # --- Panel 1: UMI counts distribution ---
    ax = axes[0]
    ax.hist(adata.obs["total_counts"], bins=100, color="steelblue", edgecolor="none")
    ax.axvline(thresholds["counts_lo"], color="red", ls="--", lw=1.5,
               label=f'low = {thresholds["counts_lo"]:.0f}')
    ax.axvline(thresholds["counts_hi"], color="red", ls="-", lw=1.5,
               label=f'high = {thresholds["counts_hi"]:.0f}')
    ax.set_xlabel("Total UMI counts")
    ax.set_ylabel("Cells")
    ax.set_title("UMI Counts Distribution")
    ax.legend(fontsize=8)

    # --- Panel 2: Gene counts distribution ---
    ax = axes[1]
    ax.hist(adata.obs["n_genes_by_counts"], bins=100, color="darkorange", edgecolor="none")
    ax.axvline(thresholds["genes_lo"], color="red", ls="--", lw=1.5,
               label=f'low = {thresholds["genes_lo"]:.0f}')
    ax.axvline(thresholds["genes_hi"], color="red", ls="-", lw=1.5,
               label=f'high = {thresholds["genes_hi"]:.0f}')
    ax.set_xlabel("Genes detected")
    ax.set_ylabel("Cells")
    ax.set_title("Gene Counts Distribution")
    ax.legend(fontsize=8)

    # --- Panel 3: MT% distribution ---
    ax = axes[2]
    ax.hist(adata.obs["pct_counts_mt"], bins=100, color="seagreen", edgecolor="none")
    ax.axvline(thresholds["mt_hi"], color="red", ls="-", lw=1.5,
               label=f'cutoff = {thresholds["mt_hi"]:.2f}%')
    ax.set_xlabel("Mitochondrial %")
    ax.set_ylabel("Cells")
    ax.set_title("MT% Distribution")
    ax.legend(fontsize=8)

    # --- Panel 4: Scatter total_counts vs n_genes, colored by MT% ---
    ax = axes[3]
    sc_plot = ax.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mt"],
        cmap="RdYlGn_r",
        s=1,
        alpha=0.5,
        rasterized=True,
    )
    fig.colorbar(sc_plot, ax=ax, label="MT %")
    ax.axvline(thresholds["counts_lo"], color="red", ls="--", lw=1, alpha=0.7)
    ax.axvline(thresholds["counts_hi"], color="red", ls="-", lw=1, alpha=0.7)
    ax.axhline(thresholds["genes_lo"], color="red", ls="--", lw=1, alpha=0.7)
    ax.axhline(thresholds["genes_hi"], color="red", ls="-", lw=1, alpha=0.7)
    ax.set_xlabel("Total UMI counts")
    ax.set_ylabel("Genes detected")
    ax.set_title("UMI vs Genes (color = MT%)")

    # --- Panel 5: Scatter UMI vs Genes, MT-passing cells only ---
    mt_pass = adata.obs["pct_counts_mt"] <= thresholds["mt_hi"]
    ax = axes[4]
    sc_plot2 = ax.scatter(
        adata.obs["total_counts"][mt_pass],
        adata.obs["n_genes_by_counts"][mt_pass],
        c=adata.obs["pct_counts_mt"][mt_pass],
        cmap="RdYlGn_r",
        s=1,
        alpha=0.5,
        rasterized=True,
    )
    fig.colorbar(sc_plot2, ax=ax, label="MT %")
    ax.axvline(thresholds["counts_lo"], color="red", ls="--", lw=1, alpha=0.7)
    ax.axvline(thresholds["counts_hi"], color="red", ls="-", lw=1, alpha=0.7)
    ax.axhline(thresholds["genes_lo"], color="red", ls="--", lw=1, alpha=0.7)
    ax.axhline(thresholds["genes_hi"], color="red", ls="-", lw=1, alpha=0.7)
    ax.set_xlabel("Total UMI counts")
    ax.set_ylabel("Genes detected")
    n_mt_pass = int(mt_pass.sum())
    ax.set_title(f"UMI vs Genes – MT-passing only ({n_mt_pass} cells)")

    # --- Panel 6 (optional): Doublet score distribution ---
    if doublet_scores is not None:
        ax = axes[5]
        ax.hist(doublet_scores, bins=100, color="mediumpurple", edgecolor="none")
        if doublet_threshold is not None:
            ax.axvline(doublet_threshold, color="red", ls="-", lw=1.5,
                       label=f"threshold = {doublet_threshold:.3f}")
        ax.set_xlabel("Scrublet doublet score")
        ax.set_ylabel("Cells")
        ax.set_title("Doublet Score Distribution")
        ax.legend(fontsize=8)

    # Hide any unused axes
    for i in range(n_panels, len(axes)):
        axes[i].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nQC plots saved to {output_path}")


def filter_cells(adata: ad.AnnData, thresholds: dict) -> ad.AnnData:
    """Apply QC filters and return the filtered AnnData."""
    n_before = adata.n_obs

    keep = (
        (adata.obs["total_counts"] >= thresholds["counts_lo"])
        & (adata.obs["total_counts"] <= thresholds["counts_hi"])
        & (adata.obs["n_genes_by_counts"] >= thresholds["genes_lo"])
        & (adata.obs["n_genes_by_counts"] <= thresholds["genes_hi"])
        & (adata.obs["pct_counts_mt"] <= thresholds["mt_hi"])
    )
    adata = adata[keep].copy()

    n_after = adata.n_obs
    print(f"\n--- QC Filtering ---")
    print(f"  Removed {n_before - n_after} / {n_before} cells "
          f"({(n_before - n_after) / n_before * 100:.1f}%)")
    print(f"  Remaining: {n_after} cells")
    return adata


def run_scrublet(adata: ad.AnnData, expected_doublet_rate: float) -> tuple[np.ndarray, float]:
    """Run Scrublet doublet detection. Returns (scores, threshold)."""
    print(f"\nRunning Scrublet (expected_doublet_rate={expected_doublet_rate}) ...")
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3)

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    n_doublets = predicted_doublets.sum()
    threshold = scrub.threshold_
    print(f"  Scrublet threshold: {threshold:.4f}")
    print(f"  Predicted doublets: {n_doublets} / {adata.n_obs} "
          f"({n_doublets / adata.n_obs * 100:.1f}%)")

    return doublet_scores, threshold


def remove_doublets(adata: ad.AnnData) -> ad.AnnData:
    """Remove cells flagged as doublets by Scrublet."""
    n_before = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    n_after = adata.n_obs
    print(f"  Removed {n_before - n_after} doublets. Remaining: {n_after} cells")
    return adata


def write_stats(adata: ad.AnnData, thresholds: dict,
                doublet_threshold: float | None, output_path: str) -> None:
    """Write a QC statistics file showing cell counts per bracket."""
    n_total = adata.n_obs
    counts = adata.obs["total_counts"].values
    genes = adata.obs["n_genes_by_counts"].values
    mt = adata.obs["pct_counts_mt"].values

    lines = []
    lines.append(f"QC Statistics Report")
    lines.append(f"Total cells (after removing zero-count): {n_total}")
    lines.append("")

    # UMI counts
    lo, hi = thresholds["counts_lo"], thresholds["counts_hi"]
    n_below = int((counts < lo).sum())
    n_pass = int(((counts >= lo) & (counts <= hi)).sum())
    n_above = int((counts > hi).sum())
    lines.append(f"--- UMI Counts (thresholds: {lo:.0f} – {hi:.0f}) ---")
    lines.append(f"  Below {lo:.0f}:          {n_below:>7} / {n_total} ({n_below/n_total*100:.1f}%)")
    lines.append(f"  {lo:.0f} – {hi:.0f} (pass): {n_pass:>7} / {n_total} ({n_pass/n_total*100:.1f}%)")
    lines.append(f"  Above {hi:.0f}:         {n_above:>7} / {n_total} ({n_above/n_total*100:.1f}%)")
    lines.append("")

    # Gene counts
    lo, hi = thresholds["genes_lo"], thresholds["genes_hi"]
    n_below = int((genes < lo).sum())
    n_pass = int(((genes >= lo) & (genes <= hi)).sum())
    n_above = int((genes > hi).sum())
    lines.append(f"--- Gene Counts (thresholds: {lo:.0f} – {hi:.0f}) ---")
    lines.append(f"  Below {lo:.0f}:          {n_below:>7} / {n_total} ({n_below/n_total*100:.1f}%)")
    lines.append(f"  {lo:.0f} – {hi:.0f} (pass): {n_pass:>7} / {n_total} ({n_pass/n_total*100:.1f}%)")
    lines.append(f"  Above {hi:.0f}:         {n_above:>7} / {n_total} ({n_above/n_total*100:.1f}%)")
    lines.append("")

    # MT%
    mt_hi = thresholds["mt_hi"]
    n_pass_mt = int((mt <= mt_hi).sum())
    n_fail_mt = int((mt > mt_hi).sum())
    lines.append(f"--- MT% (threshold: ≤ {mt_hi:.2f}%) ---")
    lines.append(f"  ≤ {mt_hi:.2f}% (pass):  {n_pass_mt:>7} / {n_total} ({n_pass_mt/n_total*100:.1f}%)")
    lines.append(f"  > {mt_hi:.2f}% (fail):  {n_fail_mt:>7} / {n_total} ({n_fail_mt/n_total*100:.1f}%)")
    lines.append("")

    # Doublets
    if doublet_threshold is not None and "doublet_score" in adata.obs.columns:
        scores = adata.obs["doublet_score"].values
        n_singlet = int((scores <= doublet_threshold).sum())
        n_doublet = int((scores > doublet_threshold).sum())
        lines.append(f"--- Doublets (threshold: {doublet_threshold:.4f}) ---")
        lines.append(f"  Singlet (pass):    {n_singlet:>7} / {n_total} ({n_singlet/n_total*100:.1f}%)")
        lines.append(f"  Doublet (fail):    {n_doublet:>7} / {n_total} ({n_doublet/n_total*100:.1f}%)")
    else:
        lines.append("--- Doublets ---")
        lines.append("  Skipped")

    report = "\n".join(lines) + "\n"
    with open(output_path, "w") as f:
        f.write(report)
    print(f"QC statistics saved to {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Single-cell preprocessing: QC, filtering, doublet removal.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("--input", required=True, help="Path to raw .h5ad file")
    p.add_argument("--output", required=True, help="Path for the cleaned .h5ad output")
    p.add_argument("--plot-output", default="qc_plots.png",
                   help="Path for the QC distribution plot (png or pdf)")
    p.add_argument("--stats-output", default="qc_stats.txt",
                   help="Path for the QC statistics text file")

    # Hard-floor filters
    p.add_argument("--min-genes", type=int, default=0, #But should be 200
                   help="Hard minimum: remove cells with fewer genes than this")
    p.add_argument("--min-counts", type=int, default=0, #But should be 500
                   help="Hard minimum: remove cells with fewer UMI counts than this")

    # MAD-based outlier detection
    p.add_argument("--n-mads", type=float, default=3.0,
                   help="Number of MADs for outlier detection on UMI/gene counts and MT%%")
    p.add_argument("--mt-cutoff", type=float, default=None,
                   help="Absolute MT%% cutoff. If not set, MAD-based detection is used")

    # Doublet removal
    p.add_argument("--expected-doublet-rate", type=float, default=0.06,
                   help="Expected doublet rate for Scrublet")
    p.add_argument("--skip-doublets", action="store_true",
                   help="Skip doublet detection entirely")

    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    # 1. Load
    adata = load_adata(args.input)

    # 2. QC metrics
    print("\nComputing QC metrics ...")
    compute_qc_metrics(adata)

    # 3. Cell cycle scoring
    print("\nScoring cell cycle ...")
    score_cell_cycle(adata)

    # 4. Thresholds
    thresholds = compute_thresholds(adata, args)

    # 5. Doublet scoring (before filtering, on full data)
    doublet_scores = None
    doublet_threshold = None
    if not args.skip_doublets:
        doublet_scores, doublet_threshold = run_scrublet(adata, args.expected_doublet_rate)

    # 6. QC plots (drawn on the unfiltered data so the full distribution is visible)
    plot_qc(adata, thresholds, doublet_scores, doublet_threshold, args.plot_output)

    # 6b. QC statistics report
    write_stats(adata, thresholds, doublet_threshold, args.stats_output)

    # 7. Filter
    adata = filter_cells(adata, thresholds)

    # 8. Remove doublets
    if not args.skip_doublets:
        adata = remove_doublets(adata)

    # 9. Save
    print(f"\nSaving cleaned object to {args.output} ...")
    adata.write_h5ad(args.output)
    print(f"  Final shape: {adata.n_obs} cells x {adata.n_vars} genes")
    print("Done.")


if __name__ == "__main__":
    main()
