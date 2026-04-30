#!/usr/bin/env python
"""
cycif_to_h5ad.py
Convert CyCIF CSV exports from cycif_export_csvs.m into an AnnData .h5ad
file compatible with cellxgene.

Usage (from command line):
    python cycif_to_h5ad.py <csv_dir> <output.h5ad>

Or called from MATLAB via cycif_to_h5ad.m.

AnnData structure produced:
    adata.X          [nCells x nMarkers]  log2 marker intensities (float32)
    adata.obs        [nCells x nMeta]     condition/treatment metadata
    adata.var        [nMarkers x 1]       marker names index
    adata.obsm       dict                 pre-computed embeddings:
                                            'X_tsne', 'X_pca', etc.
                                            (key = 'X_' + embedding name)
"""

import sys
import os
import glob
import numpy as np
import pandas as pd
import anndata as ad

def main(csv_dir, out_path):
    print(f"[cycif_to_h5ad] Reading from: {csv_dir}")

    # --- X matrix -----------------------------------------------------------
    X_path = os.path.join(csv_dir, "X.csv")
    assert os.path.exists(X_path), f"X.csv not found in {csv_dir}"
    X_df = pd.read_csv(X_path)
    X    = X_df.values.astype(np.float32)
    print(f"[cycif_to_h5ad] X: {X.shape[0]} cells x {X.shape[1]} markers")

    # --- var (marker names) -------------------------------------------------
    var_path = os.path.join(csv_dir, "var.csv")
    var_df   = pd.read_csv(var_path, index_col=0)
    var_df.index = X_df.columns.tolist()   # marker names as index

    # --- obs (cell metadata) ------------------------------------------------
    obs_path = os.path.join(csv_dir, "obs.csv")
    obs_df   = pd.read_csv(obs_path, dtype=str)

    # Convert numeric-looking columns back to numeric where appropriate
    for col in ['xy', 'cellid']:
        if col in obs_df.columns:
            obs_df[col] = pd.to_numeric(obs_df[col], errors='coerce')

    # cellxgene requires string index on obs
    obs_df.index = [f"cell_{i}" for i in range(len(obs_df))]

    # Force all object/mixed columns to explicit str first —
    # h5py cannot implicitly convert non-string objects (e.g. NaN in str cols)
    for col in obs_df.columns:
        if obs_df[col].dtype == object:
            obs_df[col] = obs_df[col].fillna('').astype(str)

    # Categorical columns improve cellxgene filter performance
    cat_cols = ['tx1_label', 'ptx_pool_label', 'ptx_label',
                'ptx_name', 'tx1_name', 'tx2_name', 'tx3_name',
                'cluster', 'condition', 'cell_type']
    for col in cat_cols:
        if col in obs_df.columns:
            obs_df[col] = obs_df[col].astype('category')

    # --- Build AnnData ------------------------------------------------------
    adata      = ad.AnnData(X=X, obs=obs_df, var=var_df)
    adata.var_names = X_df.columns.tolist()

    # Store raw so cellxgene's 'Use Raw' checkbox shows un-normalized values
    adata.raw = adata

    # --- Embeddings ---------------------------------------------------------
    emb_files = sorted(glob.glob(os.path.join(csv_dir, "emb_*.csv")))
    for ef in emb_files:
        name     = os.path.basename(ef)[4:-4]   # strip 'emb_' and '.csv'
        emb_df   = pd.read_csv(ef)
        emb_key  = f"X_{name}"
        adata.obsm[emb_key] = emb_df.values.astype(np.float32)
        print(f"[cycif_to_h5ad] Added obsm['{emb_key}'] shape={adata.obsm[emb_key].shape}")

    # --- Diagnostics before write -------------------------------------------
    n_nan  = int(np.sum(np.isnan(X)))
    n_zero = int(np.sum(X == 0))
    print(f"[cycif_to_h5ad] Pre-write  — NaN: {n_nan} ({100*n_nan/X.size:.1f}%)  "
          f"Zero: {n_zero}  Range: {np.nanmin(X):.2f} to {np.nanmax(X):.2f}")

    # --- Write .h5ad --------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    adata.write_h5ad(out_path)

    # --- Verify NaN survived round-trip -------------------------------------
    adata2  = ad.read_h5ad(out_path)
    n_nan2  = int(np.sum(np.isnan(adata2.X)))
    n_zero2 = int(np.sum(adata2.X == 0))
    print(f"[cycif_to_h5ad] Post-write — NaN: {n_nan2} ({100*n_nan2/adata2.X.size:.1f}%)  "
          f"Zero: {n_zero2}  Range: {np.nanmin(adata2.X):.2f} to {np.nanmax(adata2.X):.2f}")

    print(f"[cycif_to_h5ad] Saved: {out_path}")
    print(f"  {adata.n_obs} cells | {adata.n_vars} markers | "
          f"{len(adata.obsm)} embeddings")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python cycif_to_h5ad.py <csv_dir> <output.h5ad>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])