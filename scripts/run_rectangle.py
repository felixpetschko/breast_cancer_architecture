#!/usr/bin/env python3
"""Run Rectangle spatial deconvolution for one slide bulk matrix."""

from __future__ import annotations

import argparse
import inspect
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/rectangle.yaml")
    parser.add_argument("--reference-h5ad", default=None)
    parser.add_argument("--bulks-csv", required=True)
    parser.add_argument("--cell-type-col", default=None)
    parser.add_argument("--selected-cell-types", nargs="*", default=None)
    parser.add_argument("--correct-mrna-bias", dest="correct_mrna_bias", action="store_true")
    parser.add_argument("--no-correct-mrna-bias", dest="correct_mrna_bias", action="store_false")
    parser.set_defaults(correct_mrna_bias=None)
    parser.add_argument("--n-cpus", type=int, default=None)
    parser.add_argument("--gene-symbol-policy", choices=["none", "upper", "lower"], default=None)
    parser.add_argument("--min-genes-overlap", type=int, default=None)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--diagnostics-json", required=True)
    return parser.parse_args()


def load_config(path: str) -> dict[str, Any]:
    cfg_path = Path(path)
    if not cfg_path.exists():
        return {}
    with cfg_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def choose(cli_value: Any, cfg_value: Any, default: Any = None) -> Any:
    if cli_value is not None:
        return cli_value
    if cfg_value is not None:
        return cfg_value
    return default


def harmonize_symbols(values: list[str], policy: str) -> list[str]:
    if policy == "upper":
        return [x.upper() for x in values]
    if policy == "lower":
        return [x.lower() for x in values]
    return values


def subset_or_map_cell_types(
    ref: Any,
    cell_type_col: str,
    selected_cell_types: list[str] | None,
) -> Any:
    if cell_type_col not in ref.obs.columns:
        raise ValueError(f"Cell type column '{cell_type_col}' missing in reference obs")
    ref.obs[cell_type_col] = ref.obs[cell_type_col].astype(str)
    if selected_cell_types:
        keep = set(selected_cell_types)
        mask = ref.obs[cell_type_col].isin(keep).values
        ref = ref[mask].copy()
    return ref


def kwargs_from_signature(
    fn: Any,
    aliases: dict[str, list[str]],
    values: dict[str, Any],
) -> dict[str, Any]:
    params = set(inspect.signature(fn).parameters)
    out: dict[str, Any] = {}
    for canonical, candidate_names in aliases.items():
        if canonical not in values:
            continue
        value = values[canonical]
        if value is None:
            continue
        for name in candidate_names:
            if name in params:
                out[name] = value
                break
    return out


def run_build_signatures(rectangle: Any, payload: dict[str, Any]) -> Any:
    fn = rectangle.pp.build_rectangle_signatures
    aliases = {
        "reference": ["sc_adata", "adata_sc", "reference", "adata"],
        "bulks": ["bulks", "bulk", "bulk_df"],
        "cell_type_col": ["cell_type_col", "cell_type_column", "celltype_col", "labels_col"],
        "selected_cell_types": ["selected_cell_types", "cell_types"],
        "optimize_cutoffs": ["optimize_cutoffs"],
        "p": ["p"],
        "lfc": ["lfc"],
    }
    kwargs = kwargs_from_signature(fn, aliases, payload)
    return fn(**kwargs)


def run_deconvolution(rectangle: Any, payload: dict[str, Any]) -> Any:
    fn = rectangle.tl.deconvolution
    aliases = {
        "signatures": ["signature_result", "signatures", "signature", "sig"],
        "bulks": ["bulks", "bulk", "bulk_df"],
        "correct_mrna_bias": ["correct_mrna_bias"],
        "n_cpus": ["n_cpus", "n_jobs"],
    }
    kwargs = kwargs_from_signature(fn, aliases, payload)
    return fn(**kwargs)


def extract_fraction_table(result: Any) -> pd.DataFrame:
    if isinstance(result, pd.DataFrame):
        return result

    if isinstance(result, (list, tuple)):
        for item in result:
            try:
                return extract_fraction_table(item)
            except TypeError:
                continue

    if isinstance(result, dict):
        for key in ["fractions", "cell_fractions", "proportions", "deconv", "weights"]:
            value = result.get(key)
            if value is not None:
                try:
                    return extract_fraction_table(value)
                except TypeError:
                    continue

    for attr in ["fractions", "cell_fractions", "proportions", "deconv", "weights", "result"]:
        if hasattr(result, attr):
            value = getattr(result, attr)
            if value is not None:
                try:
                    return extract_fraction_table(value)
                except TypeError:
                    continue

    if hasattr(result, "to_pandas"):
        try:
            out = result.to_pandas()
            if isinstance(out, pd.DataFrame):
                return out
        except Exception:
            pass

    if hasattr(result, "shape") and hasattr(result, "__array__"):
        arr = result.__array__()
        if arr.ndim == 2:
            return pd.DataFrame(arr)

    raise TypeError(
        "Could not extract fractions DataFrame from Rectangle result. "
        f"Update extract_fraction_table for your rectanglepy version. Got type: {type(result)}"
    )


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    reference_path = choose(args.reference_h5ad, cfg.get("reference_path"))
    cell_type_col = choose(args.cell_type_col, cfg.get("cell_type_column"), "cell_type")
    selected_cell_types = choose(args.selected_cell_types, cfg.get("selected_cell_types"), None)
    if selected_cell_types == []:
        selected_cell_types = None
    correct_mrna_bias = choose(args.correct_mrna_bias, cfg.get("correct_mrna_bias"), True)
    n_cpus = choose(args.n_cpus, cfg.get("n_cpus"), None)
    gene_symbol_policy = choose(args.gene_symbol_policy, cfg.get("gene_symbol_policy"), "upper")
    min_genes_overlap = int(choose(args.min_genes_overlap, cfg.get("min_genes_overlap"), 500))

    if not reference_path:
        raise ValueError("Missing reference path. Set --reference-h5ad or config.reference_path")

    try:
        import anndata as ad
        import rectanglepy as rectangle
    except ImportError as exc:
        raise ImportError(
            "Missing Python dependencies. Create/activate environment from environment/rectangle.yml."
        ) from exc

    ref = ad.read_h5ad(reference_path)
    bulks = pd.read_csv(args.bulks_csv)
    if "spot_id" not in bulks.columns:
        raise ValueError("Bulk CSV must contain a 'spot_id' column")
    bulks = bulks.set_index("spot_id")

    ref.var_names = pd.Index(harmonize_symbols(ref.var_names.astype(str).tolist(), gene_symbol_policy))
    ref.var_names_make_unique()
    bulks.columns = harmonize_symbols([str(c) for c in bulks.columns], gene_symbol_policy)

    ref = subset_or_map_cell_types(ref, cell_type_col, selected_cell_types)

    overlap = sorted(set(ref.var_names).intersection(set(bulks.columns)))
    if len(overlap) < min_genes_overlap:
        raise ValueError(
            f"Gene overlap too low ({len(overlap)} < {min_genes_overlap}). "
            "Check gene symbol harmonization and reference choice."
        )

    ref = ref[:, overlap].copy()
    bulks = bulks.loc[:, overlap]

    signature_result = run_build_signatures(
        rectangle,
        {
            "reference": ref,
            "bulks": bulks,
            "cell_type_col": cell_type_col,
            "selected_cell_types": selected_cell_types,
            "optimize_cutoffs": True,
            "p": 0.01,
            "lfc": 0.5,
        },
    )

    deconv_result = run_deconvolution(
        rectangle,
        {
            "signatures": signature_result,
            "bulks": bulks,
            "correct_mrna_bias": bool(correct_mrna_bias),
            "n_cpus": n_cpus,
        },
    )

    fractions = extract_fraction_table(deconv_result).copy()
    if "spot_id" not in fractions.columns:
        fractions.index.name = "spot_id"
        fractions = fractions.reset_index()

    output_path = Path(args.out_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fractions.to_csv(output_path, index=False)

    diagnostics = {
        "reference_path": str(reference_path),
        "cell_type_col": cell_type_col,
        "selected_cell_types": selected_cell_types,
        "correct_mrna_bias": bool(correct_mrna_bias),
        "n_cpus": n_cpus,
        "gene_symbol_policy": gene_symbol_policy,
        "n_reference_cells": int(ref.n_obs),
        "n_spots": int(bulks.shape[0]),
        "n_genes_overlap": int(len(overlap)),
        "min_genes_overlap": int(min_genes_overlap),
        "fraction_columns": [c for c in fractions.columns if c != "spot_id"],
    }
    diagnostics_path = Path(args.diagnostics_json)
    diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
    diagnostics_path.write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")

    print(f"Wrote Rectangle fractions: {output_path}")
    print(f"Wrote diagnostics: {diagnostics_path}")


if __name__ == "__main__":
    main()
