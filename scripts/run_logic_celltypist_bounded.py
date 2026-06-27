#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODULE_NAME = "celltypist_bounded"
DEFAULT_MODEL_DIR = Path.home() / ".celltypist" / "data" / "models"
OUTPUT_TABLE = "celltypist_bounded_annotation.tsv"
OUTPUT_FIGURE = "umap_celltypist_bounded.png"
MODEL_STATUS_TABLE = "celltypist_model_status.tsv"


class InputPaths(NamedTuple):
    counts: Path
    genes: Path
    barcodes: Path
    umap: Path
    consensus: Path


class InputSummary(NamedTuple):
    matrix_genes: int
    matrix_cells: int
    matrix_entries: int
    genes_rows: int
    barcodes_rows: int
    umap_rows: int
    consensus_rows: int


class ModelResolution(NamedTuple):
    path: Path | None
    reason: str
    local_models: list[Path]


def append_status(out_dir: Path, module: str, status: str, detail: str = "") -> None:
    status_file = out_dir / "tables" / "module_status.tsv"
    status_file.parent.mkdir(parents=True, exist_ok=True)
    if not status_file.exists():
        status_file.write_text("module\tstatus\tdetail\n")
    clean_detail = " ".join(str(detail).replace("\t", " ").splitlines())
    with status_file.open("a") as handle:
        handle.write(f"{module}\t{status}\t{clean_detail}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run CellTypist annotation only with a bounded local model. "
            "The script never asks CellTypist to list or download remote models."
        )
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Phase output directory, for example analysis/seurat_v5_logic_run.",
    )
    parser.add_argument(
        "--model",
        default="",
        help=(
            "Optional local .pkl model path, or a local model filename/name already present "
            "under ~/.celltypist/data/models. Remote model names are not downloaded."
        ),
    )
    parser.add_argument(
        "--model-dir",
        default=DEFAULT_MODEL_DIR,
        type=Path,
        help="Local CellTypist model cache to inspect for .pkl files.",
    )
    parser.add_argument(
        "--mode",
        default="best match",
        choices=("best match", "prob match"),
        help="CellTypist prediction mode.",
    )
    parser.add_argument(
        "--p-thres",
        default=0.5,
        type=float,
        help="Probability threshold used by CellTypist in prob match mode.",
    )
    return parser.parse_args()


def phase_input_paths(out_dir: Path) -> InputPaths:
    export_dir = out_dir / "objects" / "export"
    tables_dir = out_dir / "tables"
    consensus_umap = tables_dir / "annotated_consensus_umap.tsv"
    return InputPaths(
        counts=export_dir / "counts.mtx.gz",
        genes=export_dir / "genes.tsv",
        barcodes=export_dir / "barcodes.tsv",
        umap=consensus_umap if consensus_umap.exists() else export_dir / "umap.tsv",
        consensus=tables_dir / "consensus_annotation.tsv",
    )


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing required {label}: {path}")


def read_matrix_market_shape(path: Path) -> tuple[int, int, int]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("%"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                raise ValueError(f"Invalid MatrixMarket size line in {path}: {stripped}")
            return int(parts[0]), int(parts[1]), int(parts[2])
    raise ValueError(f"MatrixMarket file has no size line: {path}")


def count_lines(path: Path) -> int:
    with path.open("rt") as handle:
        return sum(1 for _ in handle)


def inspect_inputs(paths: InputPaths) -> InputSummary:
    require_file(paths.counts, "count matrix")
    require_file(paths.genes, "gene table")
    require_file(paths.barcodes, "barcode table")
    require_file(paths.umap, "UMAP table")
    require_file(paths.consensus, "consensus annotation table")

    matrix_genes, matrix_cells, matrix_entries = read_matrix_market_shape(paths.counts)
    genes_rows = count_lines(paths.genes)
    barcodes_rows = count_lines(paths.barcodes)
    umap_rows = max(count_lines(paths.umap) - 1, 0)
    consensus_rows = max(count_lines(paths.consensus) - 1, 0)
    return InputSummary(
        matrix_genes=matrix_genes,
        matrix_cells=matrix_cells,
        matrix_entries=matrix_entries,
        genes_rows=genes_rows,
        barcodes_rows=barcodes_rows,
        umap_rows=umap_rows,
        consensus_rows=consensus_rows,
    )


def list_local_models(model_dir: Path) -> list[Path]:
    model_dir = model_dir.expanduser()
    if not model_dir.is_dir():
        return []
    return sorted(path for path in model_dir.glob("*.pkl") if path.is_file())


def is_model_file(path: Path) -> bool:
    return path.is_file() and path.suffix == ".pkl"


def resolve_requested_model(requested: str, model_dir: Path, local_models: list[Path]) -> ModelResolution:
    model_dir = model_dir.expanduser()
    requested_path = Path(requested).expanduser()

    if is_model_file(requested_path):
        return ModelResolution(requested_path, f"using requested model path {requested_path}", local_models)

    if not requested_path.is_absolute():
        cwd_path = (Path.cwd() / requested_path).resolve()
        if is_model_file(cwd_path):
            return ModelResolution(cwd_path, f"using requested model path {cwd_path}", local_models)

    request_has_path = "/" in requested or "\\" in requested or requested.startswith(".") or requested.startswith("~")
    if not request_has_path:
        wanted_names = {requested}
        if not requested.endswith(".pkl"):
            wanted_names.add(f"{requested}.pkl")
        for model_path in local_models:
            if model_path.name in wanted_names or model_path.stem == requested:
                return ModelResolution(model_path, f"using local cached model {model_path.name}", local_models)

    return ModelResolution(
        None,
        f"requested model was not found as a local .pkl file: {requested}",
        local_models,
    )


def resolve_model(model_arg: str, model_dir: Path) -> ModelResolution:
    local_models = list_local_models(model_dir)
    if model_arg:
        return resolve_requested_model(model_arg, model_dir, local_models)

    if not local_models:
        return ModelResolution(
            None,
            f"no local .pkl models found under {model_dir.expanduser()}",
            local_models,
        )

    preferred = next((path for path in local_models if path.name == "Immune_All_Low.pkl"), None)
    selected = preferred if preferred is not None else local_models[0]
    return ModelResolution(selected, f"using local cached model {selected.name}", local_models)


def load_celltypist_model(model_path: Path):
    from celltypist import models

    return models.Model.load(str(model_path.expanduser().resolve()))


def read_genes(path: Path) -> list[str]:
    genes = pd.read_csv(path, sep="\t", header=None, usecols=[0], dtype=str)[0].fillna("").astype(str).tolist()
    if not genes:
        raise ValueError(f"No genes found in {path}")
    return genes


def read_barcodes(path: Path) -> list[str]:
    barcodes = pd.read_csv(path, sep="\t", header=None, usecols=[0], dtype=str)[0].fillna("").astype(str).tolist()
    if not barcodes:
        raise ValueError(f"No barcodes found in {path}")
    if len(set(barcodes)) != len(barcodes):
        raise ValueError("barcodes.tsv contains duplicated barcodes")
    return barcodes


def read_umap(path: Path) -> pd.DataFrame:
    umap = pd.read_csv(path, sep="\t")
    if "barcode" not in umap.columns:
        raise ValueError(f"UMAP table must include a barcode column: {path}")
    if umap["barcode"].duplicated().any():
        raise ValueError("UMAP table contains duplicated barcodes")
    return umap


def read_consensus(path: Path) -> pd.DataFrame:
    consensus = pd.read_csv(path, sep="\t")
    if "barcode" not in consensus.columns:
        raise ValueError(f"Consensus annotation table must include a barcode column: {path}")
    if consensus["barcode"].duplicated().any():
        raise ValueError("Consensus annotation table contains duplicated barcodes")
    return consensus


def orient_counts_cell_by_gene(matrix, genes: list[str], barcodes: list[str]) -> sparse.csr_matrix:
    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    matrix = matrix.tocsr()

    if matrix.shape == (len(genes), len(barcodes)):
        return matrix.T.tocsr()
    if matrix.shape == (len(barcodes), len(genes)):
        return matrix.tocsr()

    raise ValueError(
        "Counts matrix shape does not match genes/barcodes: "
        f"matrix={matrix.shape}, genes={len(genes)}, barcodes={len(barcodes)}"
    )


def align_metadata_to_consensus(
    matrix: sparse.csr_matrix,
    barcodes: list[str],
    umap: pd.DataFrame,
    consensus: pd.DataFrame,
) -> tuple[sparse.csr_matrix, pd.DataFrame, pd.DataFrame]:
    barcode_index = pd.Series(np.arange(len(barcodes)), index=pd.Index(barcodes, name="barcode"))

    consensus_obs = consensus.copy()
    consensus_obs["barcode"] = consensus_obs["barcode"].astype(str)
    consensus_obs["_matrix_row"] = consensus_obs["barcode"].map(barcode_index)
    missing = int(consensus_obs["_matrix_row"].isna().sum())
    if missing:
        raise ValueError(f"{missing} consensus barcodes are missing from barcodes.tsv")

    consensus_obs["_matrix_row"] = consensus_obs["_matrix_row"].astype(int)
    matrix = matrix[consensus_obs["_matrix_row"].to_numpy(), :].tocsr()
    obs = consensus_obs.drop(columns=["_matrix_row"]).set_index("barcode", drop=False)

    umap_obs = umap.copy()
    umap_obs["barcode"] = umap_obs["barcode"].astype(str)
    umap_obs = pd.DataFrame({"barcode": obs["barcode"].to_numpy()}).merge(umap_obs, on="barcode", how="left")
    return matrix, obs, umap_obs


def build_adata(paths: InputPaths):
    import scanpy as sc
    from anndata import AnnData

    genes = read_genes(paths.genes)
    barcodes = read_barcodes(paths.barcodes)
    umap = read_umap(paths.umap)
    consensus = read_consensus(paths.consensus)
    matrix = orient_counts_cell_by_gene(mmread(paths.counts), genes, barcodes).astype(np.float32)

    if matrix.shape[1] != len(genes):
        raise ValueError(f"Gene count mismatch: matrix has {matrix.shape[1]} genes, genes.tsv has {len(genes)} rows")
    if matrix.shape[0] != len(barcodes):
        raise ValueError(f"Cell count mismatch: matrix has {matrix.shape[0]} cells, barcodes.tsv has {len(barcodes)} rows")

    matrix, obs, umap = align_metadata_to_consensus(matrix, barcodes, umap, consensus)

    var = pd.DataFrame(index=pd.Index(pd.Series(genes, dtype=str), name="gene"))
    adata = AnnData(X=matrix, obs=obs, var=var)

    if {"umap_1", "umap_2"}.issubset(umap.columns):
        adata.obsm["X_umap"] = umap[["umap_1", "umap_2"]].to_numpy(dtype=float)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata, umap, consensus


def run_celltypist(adata, model, mode: str, p_thres: float):
    import celltypist

    return celltypist.annotate(
        adata,
        model=model,
        mode=mode,
        p_thres=p_thres,
        majority_voting=False,
    )


def predictions_to_table(result, adata, model_path: Path) -> pd.DataFrame:
    barcodes = pd.Index(adata.obs_names.astype(str), name="barcode")
    labels = result.predicted_labels.copy()
    labels.index = labels.index.astype(str)
    labels = labels.reindex(barcodes)

    probability = result.probability_matrix.copy()
    probability.index = probability.index.astype(str)
    probability = probability.reindex(barcodes)

    output = pd.DataFrame({"barcode": barcodes.to_numpy()})
    if "predicted_labels" in labels.columns:
        output["celltypist_label"] = labels["predicted_labels"].astype(str).to_numpy()
    elif labels.shape[1] > 0:
        output["celltypist_label"] = labels.iloc[:, 0].astype(str).to_numpy()
    else:
        output["celltypist_label"] = ""

    for column in labels.columns:
        output[f"celltypist_{column}"] = labels[column].astype(str).to_numpy()

    if probability.shape[1] > 0:
        output["celltypist_conf_score"] = probability.max(axis=1).to_numpy(dtype=float)
        output["celltypist_best_probability_label"] = probability.idxmax(axis=1).astype(str).to_numpy()
    else:
        output["celltypist_conf_score"] = np.nan
        output["celltypist_best_probability_label"] = ""

    for column in ("sample_id", "cluster", "logic_marker_rule_celltype", "logic_consensus_celltype"):
        if column in adata.obs.columns:
            output[column] = adata.obs[column].astype("string").fillna("").to_numpy()

    output["celltypist_model"] = model_path.name
    return output


def write_skip_table(
    out_dir: Path,
    reason: str,
    model_arg: str,
    model_dir: Path,
    local_models: list[Path],
    summary: InputSummary,
) -> None:
    row = {
        "status": "SKIPPED_DEPENDENCY",
        "reason": reason,
        "model_requested": model_arg,
        "model_dir": str(model_dir.expanduser()),
        "local_pkl_models": len(local_models),
        "counts_genes": summary.matrix_genes,
        "counts_cells": summary.matrix_cells,
        "genes_rows": summary.genes_rows,
        "barcodes_rows": summary.barcodes_rows,
        "umap_rows": summary.umap_rows,
        "consensus_rows": summary.consensus_rows,
    }
    pd.DataFrame([row]).to_csv(out_dir / "tables" / OUTPUT_TABLE, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "module": MODULE_NAME,
                "status": "SKIPPED_DEPENDENCY",
                "model_requested": model_arg,
                "model_path": "",
                "model_dir": str(model_dir.expanduser()),
                "local_pkl_models": len(local_models),
                "reason": reason,
                "output_table": f"tables/{OUTPUT_TABLE}",
            }
        ]
    ).to_csv(out_dir / "tables" / MODEL_STATUS_TABLE, sep="\t", index=False)


def write_pass_model_status(
    out_dir: Path,
    model_arg: str,
    model_dir: Path,
    model_path: Path,
    local_models: list[Path],
    detail: str,
) -> None:
    pd.DataFrame(
        [
            {
                "module": MODULE_NAME,
                "status": "PASS",
                "model_requested": model_arg,
                "model_path": str(model_path),
                "model_dir": str(model_dir.expanduser()),
                "local_pkl_models": len(local_models),
                "reason": detail,
                "output_table": f"tables/{OUTPUT_TABLE}",
            }
        ]
    ).to_csv(out_dir / "tables" / MODEL_STATUS_TABLE, sep="\t", index=False)


def plot_umap(predictions: pd.DataFrame, umap: pd.DataFrame, out_dir: Path) -> bool:
    required = {"barcode", "umap_1", "umap_2"}
    if predictions.empty or "celltypist_label" not in predictions.columns or not required.issubset(umap.columns):
        return False

    plot_df = umap[["barcode", "umap_1", "umap_2"]].merge(
        predictions[["barcode", "celltypist_label"]],
        on="barcode",
        how="inner",
    )
    plot_df = plot_df.dropna(subset=["umap_1", "umap_2", "celltypist_label"])
    if plot_df.empty:
        return False

    label_counts = plot_df["celltypist_label"].value_counts()
    labels = label_counts.index.tolist()
    cmap = plt.get_cmap("tab20")
    color_by_label = {label: cmap(i % cmap.N) for i, label in enumerate(labels)}

    fig, ax = plt.subplots(figsize=(7, 5))
    for label in labels:
        sub = plot_df[plot_df["celltypist_label"] == label]
        legend_label = f"{label} ({len(sub)})" if len(labels) <= 15 else None
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=4,
            linewidths=0,
            alpha=0.8,
            c=[color_by_label[label]],
            label=legend_label,
        )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("CellTypist bounded annotation")
    if len(labels) <= 15:
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=3)
    fig.tight_layout()
    fig.savefig(out_dir / "figures" / OUTPUT_FIGURE, dpi=220)
    plt.close(fig)
    return True


def main() -> None:
    args = parse_args()
    out_dir = args.out
    model_dir = args.model_dir
    out_dir.joinpath("tables").mkdir(parents=True, exist_ok=True)
    out_dir.joinpath("figures").mkdir(parents=True, exist_ok=True)

    append_status(out_dir, MODULE_NAME, "START", "bounded local CellTypist model resolution")
    try:
        paths = phase_input_paths(out_dir)
        summary = inspect_inputs(paths)

        resolution = resolve_model(args.model, model_dir)
        if resolution.path is None:
            write_skip_table(out_dir, resolution.reason, args.model, model_dir, resolution.local_models, summary)
            append_status(out_dir, MODULE_NAME, "SKIPPED_DEPENDENCY", resolution.reason)
            print(f"SKIPPED_DEPENDENCY: {resolution.reason}")
            return

        try:
            model = load_celltypist_model(resolution.path)
        except Exception as exc:
            reason = f"local model could not be loaded: {type(exc).__name__}: {exc}"
            write_skip_table(out_dir, reason, args.model, model_dir, resolution.local_models, summary)
            append_status(out_dir, MODULE_NAME, "SKIPPED_DEPENDENCY", reason)
            print(f"SKIPPED_DEPENDENCY: {reason}")
            return

        adata, umap, _consensus = build_adata(paths)
        result = run_celltypist(adata, model, args.mode, args.p_thres)
        predictions = predictions_to_table(result, adata, resolution.path)
        predictions.to_csv(out_dir / "tables" / OUTPUT_TABLE, sep="\t", index=False)
        plotted = plot_umap(predictions, umap, out_dir)

        label_count = predictions["celltypist_label"].nunique(dropna=True)
        detail = (
            f"{len(predictions)} cells annotated; labels={label_count}; "
            f"model={resolution.path.name}; figure={'yes' if plotted else 'no'}"
        )
        write_pass_model_status(out_dir, args.model, model_dir, resolution.path, resolution.local_models, detail)
        append_status(out_dir, MODULE_NAME, "PASS", detail)
        print(detail)
    except Exception as exc:
        append_status(out_dir, MODULE_NAME, "FAIL", f"{type(exc).__name__}: {exc}")
        raise


if __name__ == "__main__":
    main()
