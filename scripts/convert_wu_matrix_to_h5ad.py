#!/usr/bin/env python3
"""Convert Wu et al. GEO mtx/tsv/csv files into AnnData h5ad."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", required=True, help="Directory with mtx/genes/barcodes/metadata")
    p.add_argument("--output-h5ad", required=True, help="Output h5ad path")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    base = Path(args.input_dir)
    out = Path(args.output_h5ad)

    mtx_path = base / "count_matrix_sparse.mtx"
    genes_path = base / "count_matrix_genes.tsv"
    barcodes_path = base / "count_matrix_barcodes.tsv"
    meta_path = base / "metadata.csv"

    for p in [mtx_path, genes_path, barcodes_path, meta_path]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required input file: {p}")

    print("Loading sparse matrix from", mtx_path)
    # Matrix Market is genes x cells. Convert to cells x genes for AnnData.
    X = scipy.io.mmread(mtx_path).tocsr().T.tocsr()
    print(f"Matrix loaded: shape={X.shape}, nnz={X.nnz}")

    print("Loading genes/barcodes/metadata")
    genes = pd.read_csv(genes_path, header=None, sep="\t")[0].astype(str)
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].astype(str)
    meta = pd.read_csv(meta_path, index_col=0)
    meta.index = meta.index.astype(str)

    if X.shape[0] != len(barcodes):
        raise ValueError(f"Cells mismatch: matrix={X.shape[0]}, barcodes={len(barcodes)}")
    if X.shape[1] != len(genes):
        raise ValueError(f"Genes mismatch: matrix={X.shape[1]}, genes={len(genes)}")

    obs = pd.DataFrame(index=barcodes.values)
    obs = obs.join(meta, how="left")
    missing_meta = int(obs.isna().all(axis=1).sum())

    var = pd.DataFrame(index=genes.values)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)

    print("Wrote", out)
    print("AnnData:", adata)
    print("Cells without metadata:", missing_meta)


if __name__ == "__main__":
    main()
