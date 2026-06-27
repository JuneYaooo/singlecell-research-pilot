#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MODULE_NAME = "communication_lr"
CELLTYPE_COLUMNS = (
    "logic_consensus_celltype",
    "consensus_celltype",
    "celltype",
    "annotation",
)
SCORE_COLUMNS = [
    "ligand",
    "receptor",
    "receptor_genes",
    "sender_celltype",
    "receiver_celltype",
    "sender_n_cells",
    "receiver_n_cells",
    "avg_ligand_sender",
    "avg_receptor_receiver",
    "score",
    "ligand_present",
    "receptor_genes_present",
    "missing_genes",
]


@dataclass(frozen=True)
class LRPair:
    ligand: str
    receptors: tuple[str, ...]

    @property
    def receptor_label(self) -> str:
        return "/".join(self.receptors)


CURATED_LR_PAIRS = (
    LRPair("CXCL12", ("CXCR4",)),
    LRPair("CCL5", ("CCR5",)),
    LRPair("IL7", ("IL7R",)),
    LRPair("TGFB1", ("TGFBR1", "TGFBR2")),
    LRPair("VEGFA", ("KDR", "FLT1")),
    LRPair("COL1A1", ("ITGA1", "ITGB1")),
    LRPair("SPP1", ("CD44",)),
    LRPair("CD40LG", ("CD40",)),
    LRPair("LTA", ("TNFRSF1A",)),
    LRPair("ICAM1", ("ITGAL",)),
)


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
            "Compute local curated ligand-receptor communication scores from "
            "Phase 1 exported counts and consensus annotations."
        )
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Phase 1 output directory, for example analysis/seurat_v5_logic_run.",
    )
    parser.add_argument(
        "--min-cells",
        default=20,
        type=int,
        help="Minimum cells required for a cell type to enter communication scoring.",
    )
    return parser.parse_args()


def clean_celltype(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "na"}:
        return ""
    return text


def normalize_gene(value: object) -> str:
    return str(value).strip().upper()


def curated_genes() -> set[str]:
    genes: set[str] = set()
    for pair in CURATED_LR_PAIRS:
        genes.add(pair.ligand)
        genes.update(pair.receptors)
    return genes


def read_genes(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Missing genes file: {path}")

    table = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if table.empty:
        raise ValueError(f"Genes file is empty: {path}")

    expected = curated_genes()
    best_column = 0
    best_hits = -1
    for column in table.columns:
        hits = table[column].map(normalize_gene).isin(expected).sum()
        if hits > best_hits:
            best_column = column
            best_hits = int(hits)

    genes = table[best_column].fillna("").astype(str).str.strip().tolist()
    if not genes:
        raise ValueError(f"Genes file has no gene names: {path}")
    return genes


def read_barcodes(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Missing barcode file: {path}")
    table = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if table.empty:
        raise ValueError(f"Barcode file is empty: {path}")
    barcodes = table.iloc[:, 0].fillna("").astype(str).str.strip().tolist()
    if len(set(barcodes)) != len(barcodes):
        raise ValueError("Barcode file contains duplicate barcodes")
    return barcodes


def read_annotations(path: Path) -> tuple[pd.DataFrame, str]:
    if not path.exists():
        raise FileNotFoundError(f"Missing annotation file: {path}")

    annotations = pd.read_csv(path, sep="\t")
    if annotations.empty:
        raise ValueError(f"Annotation file is empty: {path}")

    celltype_column = next((col for col in CELLTYPE_COLUMNS if col in annotations.columns), None)
    if celltype_column is None:
        raise ValueError(
            "Annotation file must include one of these cell type columns: "
            + ", ".join(CELLTYPE_COLUMNS)
        )

    annotations = annotations.copy()
    annotations["_celltype"] = annotations[celltype_column].map(clean_celltype)
    return annotations, celltype_column


def read_counts(path: Path) -> sparse.csr_matrix:
    if not path.exists():
        raise FileNotFoundError(f"Missing counts matrix: {path}")

    matrix = mmread(str(path))
    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    return matrix.tocsr()


def orient_counts(matrix: sparse.spmatrix, n_genes: int, n_barcodes: int) -> sparse.csr_matrix:
    if matrix.shape == (n_genes, n_barcodes):
        return matrix.tocsr()
    if matrix.shape == (n_barcodes, n_genes):
        return matrix.T.tocsr()
    raise ValueError(
        "Counts matrix shape does not match genes and barcodes: "
        f"matrix={matrix.shape}, genes={n_genes}, barcodes={n_barcodes}"
    )


def align_counts_to_annotations(
    counts: sparse.csr_matrix,
    annotations: pd.DataFrame,
    barcodes: list[str],
) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    barcode_index = pd.Series(np.arange(len(barcodes)), index=barcodes)
    aligned = annotations.copy()
    aligned["_matrix_index"] = aligned["barcode"].astype(str).map(barcode_index)
    missing = int(aligned["_matrix_index"].isna().sum())
    if missing:
        aligned = aligned.dropna(subset=["_matrix_index"]).copy()
    if aligned.empty:
        raise ValueError("No consensus annotation barcodes were found in the exported count matrix")
    aligned["_matrix_index"] = aligned["_matrix_index"].astype(int)
    aligned = aligned.sort_values("_matrix_index").reset_index(drop=True)
    subset = counts[:, aligned["_matrix_index"].to_numpy()].tocsr()
    aligned = aligned.drop(columns=["_matrix_index"])
    return subset, aligned


def gene_index_map(genes: list[str]) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for row_index, gene in enumerate(genes):
        key = normalize_gene(gene)
        if not key:
            continue
        index.setdefault(key, []).append(row_index)
    return index


def eligible_celltypes(annotations: pd.DataFrame, min_cells: int) -> tuple[list[str], dict[str, np.ndarray]]:
    counts = annotations.loc[annotations["_celltype"] != "", "_celltype"].value_counts()
    eligible = sorted(counts[counts >= min_cells].index.tolist())
    celltypes = annotations["_celltype"].to_numpy()
    indices = {celltype: np.flatnonzero(celltypes == celltype) for celltype in eligible}
    return eligible, indices


def average_expression(
    matrix: sparse.csr_matrix,
    gene_rows: list[int],
    cell_indices: np.ndarray,
) -> float:
    if not gene_rows or len(cell_indices) == 0:
        return 0.0
    total = matrix[gene_rows, :][:, cell_indices].sum()
    return float(total) / float(len(cell_indices))


def compute_gene_means(
    matrix: sparse.csr_matrix,
    gene_index: dict[str, list[int]],
    cell_indices_by_type: dict[str, np.ndarray],
) -> dict[str, dict[str, float]]:
    means: dict[str, dict[str, float]] = {}
    for gene in sorted(curated_genes()):
        rows = gene_index.get(gene, [])
        means[gene] = {
            celltype: average_expression(matrix, rows, indices)
            for celltype, indices in cell_indices_by_type.items()
        }
    return means


def compute_lr_scores(
    matrix: sparse.csr_matrix,
    genes: list[str],
    annotations: pd.DataFrame,
    celltypes: list[str],
    cell_indices_by_type: dict[str, np.ndarray],
) -> pd.DataFrame:
    gene_index = gene_index_map(genes)
    gene_means = compute_gene_means(matrix, gene_index, cell_indices_by_type)
    rows: list[dict[str, object]] = []

    for pair in CURATED_LR_PAIRS:
        ligand = normalize_gene(pair.ligand)
        receptors = tuple(normalize_gene(receptor) for receptor in pair.receptors)
        ligand_present = ligand in gene_index
        receptor_present = [receptor for receptor in receptors if receptor in gene_index]
        pair_genes = (ligand,) + receptors
        missing = [gene for gene in pair_genes if gene not in gene_index]

        for sender in celltypes:
            ligand_average = gene_means[ligand][sender]
            sender_n = int(len(cell_indices_by_type[sender]))
            for receiver in celltypes:
                receptor_values = [gene_means[receptor][receiver] for receptor in receptors]
                receptor_average = float(np.mean(receptor_values)) if receptor_values else 0.0
                receiver_n = int(len(cell_indices_by_type[receiver]))
                rows.append(
                    {
                        "ligand": ligand,
                        "receptor": "/".join(receptors),
                        "receptor_genes": ",".join(receptors),
                        "sender_celltype": sender,
                        "receiver_celltype": receiver,
                        "sender_n_cells": sender_n,
                        "receiver_n_cells": receiver_n,
                        "avg_ligand_sender": ligand_average,
                        "avg_receptor_receiver": receptor_average,
                        "score": ligand_average * receptor_average,
                        "ligand_present": ligand_present,
                        "receptor_genes_present": ",".join(receptor_present),
                        "missing_genes": ",".join(missing),
                    }
                )

    scores = pd.DataFrame(rows, columns=SCORE_COLUMNS)
    if scores.empty:
        return scores
    return scores.sort_values(
        ["score", "ligand", "receptor", "sender_celltype", "receiver_celltype"],
        ascending=[False, True, True, True, True],
    ).reset_index(drop=True)


def write_empty_scores(path: Path) -> None:
    pd.DataFrame(columns=SCORE_COLUMNS).to_csv(path, sep="\t", index=False)


def plot_skip_heatmap(path: Path, reason: str) -> None:
    plt.figure(figsize=(7, 3.5))
    plt.axis("off")
    plt.text(
        0.5,
        0.55,
        "Communication LR scoring skipped",
        ha="center",
        va="center",
        fontsize=13,
        weight="bold",
    )
    plt.text(0.5, 0.38, reason, ha="center", va="center", fontsize=10, wrap=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def plot_lr_heatmap(scores: pd.DataFrame, celltypes: list[str], path: Path) -> None:
    heatmap = (
        scores.groupby(["sender_celltype", "receiver_celltype"], observed=True)["score"]
        .sum()
        .unstack(fill_value=0.0)
        .reindex(index=celltypes, columns=celltypes, fill_value=0.0)
    )
    values = heatmap.to_numpy(dtype=float)

    width = max(6.0, 0.62 * len(celltypes) + 2.4)
    height = max(5.0, 0.58 * len(celltypes) + 1.8)
    fig, ax = plt.subplots(figsize=(width, height))
    image = ax.imshow(values, cmap="magma", aspect="auto")
    ax.set_xticks(np.arange(len(celltypes)))
    ax.set_yticks(np.arange(len(celltypes)))
    ax.set_xticklabels(celltypes, rotation=45, ha="right")
    ax.set_yticklabels(celltypes)
    ax.set_xlabel("Receiver cell type")
    ax.set_ylabel("Sender cell type")
    ax.set_title("Curated ligand-receptor communication score")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Sum of LR scores")

    if len(celltypes) <= 8:
        max_value = float(np.nanmax(values)) if values.size else 0.0
        threshold = max_value * 0.5
        for row in range(values.shape[0]):
            for col in range(values.shape[1]):
                value = values[row, col]
                color = "white" if value > threshold and max_value > 0 else "black"
                ax.text(col, row, f"{value:.2g}", ha="center", va="center", color=color, fontsize=8)

    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run(out_dir: Path, min_cells: int) -> tuple[str, str]:
    if min_cells < 1:
        raise ValueError("--min-cells must be at least 1")

    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    export_dir = out_dir / "objects" / "export"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    scores_path = tables_dir / "communication_lr_scores.tsv"
    heatmap_path = figures_dir / "communication_lr_heatmap.png"

    genes = read_genes(export_dir / "genes.tsv")
    barcodes = read_barcodes(export_dir / "barcodes.tsv")
    annotations, celltype_column = read_annotations(tables_dir / "consensus_annotation.tsv")
    counts = orient_counts(read_counts(export_dir / "counts.mtx.gz"), len(genes), len(barcodes))
    counts, annotations = align_counts_to_annotations(counts, annotations, barcodes)

    celltypes, cell_indices_by_type = eligible_celltypes(annotations, min_cells)
    if len(celltypes) < 2:
        available = annotations.loc[annotations["_celltype"] != "", "_celltype"].value_counts().to_dict()
        reason = (
            f"Need at least 2 cell types with >= {min_cells} cells; "
            f"found {len(celltypes)} eligible from {len(available)} annotated cell types."
        )
        write_empty_scores(scores_path)
        plot_skip_heatmap(heatmap_path, reason)
        return "SKIPPED_ASSUMPTION", reason

    scores = compute_lr_scores(counts, genes, annotations, celltypes, cell_indices_by_type)
    scores.to_csv(scores_path, sep="\t", index=False)
    plot_lr_heatmap(scores, celltypes, heatmap_path)

    nonzero = int((scores["score"] > 0).sum())
    present_genes = len(curated_genes().intersection(gene_index_map(genes)))
    detail = (
        f"{len(scores)} LR scores across {len(celltypes)} cell types using "
        f"{celltype_column}; nonzero_scores={nonzero}; curated_genes_present={present_genes}"
    )
    return "PASS", detail


def main() -> None:
    args = parse_args()
    out_dir = args.out
    append_status(out_dir, MODULE_NAME, "START", "local_curated_lr")
    try:
        status, detail = run(out_dir, args.min_cells)
    except Exception as exc:
        append_status(out_dir, MODULE_NAME, "FAIL", f"{type(exc).__name__}: {exc}")
        raise
    append_status(out_dir, MODULE_NAME, status, detail)
    print(detail)


if __name__ == "__main__":
    main()
