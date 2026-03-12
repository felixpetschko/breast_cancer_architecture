#!/usr/bin/env python3
"""Prepare Wu et al. scRNA reference for Rectangle.

This script normalizes metadata and gene symbols, optionally remaps cell labels,
and writes an AnnData h5ad file ready for Rectangle signature building.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", required=True, help="Input Wu atlas AnnData (.h5ad)")
    parser.add_argument("--output-h5ad", required=True, help="Output AnnData path")
    parser.add_argument(
        "--cell-type-col",
        required=True,
        help="obs column containing cell-type labels",
    )
    parser.add_argument(
        "--gene-symbol-policy",
        choices=["none", "upper", "lower"],
        default="upper",
        help="How to harmonize gene symbols",
    )
    parser.add_argument(
        "--label-map-csv",
        default=None,
        help="Optional CSV with columns: source_label,target_label",
    )
    parser.add_argument(
        "--selected-cell-types",
        nargs="*",
        default=None,
        help="Optional list of labels to keep after mapping",
    )
    parser.add_argument(
        "--summary-json",
        default=None,
        help="Optional path for summary diagnostics JSON",
    )
    return parser.parse_args()


def harmonize_genes(index: pd.Index, policy: str) -> pd.Index:
    genes = index.astype(str)
    if policy == "upper":
        genes = genes.str.upper()
    elif policy == "lower":
        genes = genes.str.lower()
    return pd.Index(genes)


def load_label_map(path: str | None) -> dict[str, str]:
    if path is None:
        return {}
    df = pd.read_csv(path)
    required = {"source_label", "target_label"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"label map missing columns: {sorted(missing)}")
    return dict(zip(df["source_label"].astype(str), df["target_label"].astype(str)))


def main() -> None:
    args = parse_args()

    adata = ad.read_h5ad(args.input_h5ad)

    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(f"Cell type column '{args.cell_type_col}' not found in obs")

    # Keep raw counts in X; Rectangle expects non-log-transformed count-like input.
    if np.any(np.asarray(adata.X.sum(axis=1)).ravel() <= 0):
        print("Warning: some cells have zero total counts.")

    label_map = load_label_map(args.label_map_csv)
    labels = adata.obs[args.cell_type_col].astype(str)
    if label_map:
        labels = labels.map(lambda x: label_map.get(x, x))
    adata.obs[args.cell_type_col] = labels

    genes = harmonize_genes(adata.var_names, args.gene_symbol_policy)
    adata.var_names = genes
    adata.var_names_make_unique()

    if args.selected_cell_types:
        keep = set(args.selected_cell_types)
        mask = adata.obs[args.cell_type_col].isin(keep).values
        adata = adata[mask].copy()

    out_path = Path(args.output_h5ad)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)

    summary = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "cell_type_col": args.cell_type_col,
        "n_cell_types": int(adata.obs[args.cell_type_col].nunique()),
        "cell_type_counts": adata.obs[args.cell_type_col].value_counts().to_dict(),
        "gene_symbol_policy": args.gene_symbol_policy,
    }
    if args.summary_json:
        summary_path = Path(args.summary_json)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2))

    print(f"Wrote prepared reference: {out_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
