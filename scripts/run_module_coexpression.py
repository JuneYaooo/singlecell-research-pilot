#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.io import mmread
from scipy.spatial.distance import squareform


MODULE = "coexpression"
RANDOM_SEED = 17
MAX_VARIABLE_GENES = 2000
CELLTYPE_COLUMNS = (
    "logic_consensus_celltype",
    "consensus_celltype",
    "celltype",
    "annotation",
    "logic_marker_rule_celltype",
)
SUBGROUP_COLUMNS = ("sample_id", "cluster")


@dataclass(frozen=True)
class InputPaths:
    counts: Path
    genes: Path
    barcodes: Path
    annotations: Path


@dataclass(frozen=True)
class Metacells:
    profiles: np.ndarray
    labels: list[str]
    celltypes: list[str]


def append_status(out_dir: Path, module: str, status: str, detail: str = "") -> None:
    status_file = out_dir / "tables" / "module_status.tsv"
    status_file.parent.mkdir(parents=True, exist_ok=True)
    if not status_file.exists():
        status_file.write_text("module\tstatus\tdetail\n")
    clean_detail = " ".join(str(detail).replace("\t", " ").splitlines())
    with status_file.open("a") as handle:
        handle.write(f"{module}\t{status}\t{clean_detail}\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a bounded local hdWGCNA-style co-expression replacement from "
            "core workflow exported counts and consensus annotations."
        )
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="core workflow output directory, for example analysis/workflow_run.",
    )
    parser.add_argument("--counts", type=Path, help="MatrixMarket counts path. Defaults under --out.")
    parser.add_argument("--genes", type=Path, help="genes.tsv path. Defaults under --out.")
    parser.add_argument("--barcodes", type=Path, help="barcodes.tsv path. Defaults under --out.")
    parser.add_argument(
        "--annotations",
        type=Path,
        help="consensus_annotation.tsv path. Defaults under --out/tables.",
    )
    parser.add_argument(
        "--n-variable-genes",
        default=1500,
        type=int,
        help=f"High-variance gene cap before co-expression clustering; hard-capped at {MAX_VARIABLE_GENES}.",
    )
    parser.add_argument(
        "--normalization-scale",
        default=10000.0,
        type=float,
        help="Per-cell count depth used before log1p normalization.",
    )
    parser.add_argument(
        "--min-detected-cells",
        default=5,
        type=int,
        help="Minimum cells with nonzero normalized expression for a gene to be considered variable.",
    )
    parser.add_argument(
        "--min-cells-per-metacell",
        default=20,
        type=int,
        help="Preferred minimum cells in natural metacell subgroups.",
    )
    parser.add_argument(
        "--max-cells-per-metacell",
        default=120,
        type=int,
        help="Large natural subgroups are split into deterministic chunks of about this size.",
    )
    parser.add_argument(
        "--min-module-size",
        default=25,
        type=int,
        help="Preferred minimum module size when merging small hierarchical clusters.",
    )
    parser.add_argument(
        "--max-modules",
        default=12,
        type=int,
        help="Maximum number of co-expression modules to request from hierarchical clustering.",
    )
    return parser.parse_args(argv)


def default_counts_path(export_dir: Path) -> Path:
    plain = export_dir / "counts.mtx"
    gzipped = export_dir / "counts.mtx.gz"
    return plain if plain.exists() else gzipped


def resolve_paths(args: argparse.Namespace) -> InputPaths:
    export_dir = args.out / "objects" / "export"
    return InputPaths(
        counts=args.counts or default_counts_path(export_dir),
        genes=args.genes or export_dir / "genes.tsv",
        barcodes=args.barcodes or export_dir / "barcodes.tsv",
        annotations=args.annotations or args.out / "tables" / "consensus_annotation.tsv",
    )


def require_files(paths: InputPaths) -> None:
    missing = [str(path) for path in (paths.counts, paths.genes, paths.barcodes, paths.annotations) if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input file(s): " + ", ".join(missing))


def clean_label(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "na"}:
        return ""
    return text


def read_vector(path: Path, name: str) -> list[str]:
    table = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if table.empty:
        raise ValueError(f"{name} is empty: {path}")
    values = table.iloc[:, 0].fillna("").astype(str).str.strip().tolist()
    if not values:
        raise ValueError(f"{name} has no values: {path}")
    return values


def read_barcodes(path: Path) -> list[str]:
    barcodes = read_vector(path, "barcodes.tsv")
    if len(set(barcodes)) != len(barcodes):
        raise ValueError("barcodes.tsv contains duplicate barcodes")
    return barcodes


def choose_celltype_column(annotations: pd.DataFrame) -> str:
    for column in CELLTYPE_COLUMNS:
        if column in annotations.columns:
            return column
    raise ValueError(
        "consensus_annotation.tsv missing a usable cell type column "
        f"({', '.join(CELLTYPE_COLUMNS)})"
    )


def read_annotations(path: Path) -> tuple[pd.DataFrame, str]:
    annotations = pd.read_csv(path, sep="\t", dtype=str)
    if annotations.empty:
        raise ValueError(f"Annotation file is empty: {path}")
    if "barcode" not in annotations.columns:
        raise ValueError(f"{path} must include a barcode column")

    celltype_column = choose_celltype_column(annotations)
    keep_columns = ["barcode", celltype_column]
    keep_columns.extend(column for column in SUBGROUP_COLUMNS if column in annotations.columns)

    clean = annotations[keep_columns].copy()
    clean = clean.rename(columns={celltype_column: "celltype"})
    clean["barcode"] = clean["barcode"].map(clean_label)
    clean["celltype"] = clean["celltype"].map(clean_label)
    for column in SUBGROUP_COLUMNS:
        if column in clean.columns:
            clean[column] = clean[column].map(clean_label)

    clean = clean[(clean["barcode"] != "") & (clean["celltype"] != "")].copy()
    clean = clean.drop_duplicates(subset=["barcode"], keep="first")
    if clean.empty:
        raise ValueError("No annotated singlets with usable barcode and cell type values")
    return clean, celltype_column


def read_counts(path: Path, genes: list[str], barcodes: list[str]) -> sparse.csr_matrix:
    matrix = mmread(str(path))
    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    matrix = matrix.tocsr()

    if matrix.shape == (len(genes), len(barcodes)):
        return matrix
    if matrix.shape == (len(barcodes), len(genes)):
        return matrix.T.tocsr()

    raise ValueError(
        "Counts matrix shape does not match genes/barcodes: "
        f"matrix={matrix.shape}, genes={len(genes)}, barcodes={len(barcodes)}"
    )


def align_to_annotated_singlets(
    gene_by_cell: sparse.csr_matrix,
    annotations: pd.DataFrame,
    barcodes: list[str],
) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    barcode_index = pd.Series(np.arange(len(barcodes)), index=barcodes)
    aligned = annotations.copy()
    aligned["_matrix_col"] = aligned["barcode"].map(barcode_index)
    aligned = aligned.dropna(subset=["_matrix_col"]).copy()
    if aligned.empty:
        raise ValueError("No consensus annotation barcodes were found in the exported count matrix")

    aligned["_matrix_col"] = aligned["_matrix_col"].astype(int)
    aligned = aligned.sort_values("_matrix_col").reset_index(drop=True)
    matrix = gene_by_cell[:, aligned["_matrix_col"].to_numpy()].tocsr()
    aligned = aligned.drop(columns=["_matrix_col"])
    return matrix, aligned


def normalize_log1p_gene_by_cell(gene_by_cell: sparse.csr_matrix, scale: float) -> sparse.csr_matrix:
    if scale <= 0:
        raise ValueError("--normalization-scale must be positive")

    totals = np.asarray(gene_by_cell.sum(axis=0)).ravel().astype(np.float64)
    factors = np.zeros_like(totals)
    nonzero = totals > 0
    factors[nonzero] = float(scale) / totals[nonzero]

    normalized = gene_by_cell.multiply(factors).tocsr()
    normalized.data = np.log1p(normalized.data)
    return normalized


def select_variable_genes(
    normalized: sparse.csr_matrix,
    n_variable_genes: int,
    min_detected_cells: int,
) -> np.ndarray:
    if n_variable_genes <= 0:
        raise ValueError("--n-variable-genes must be positive")
    if min_detected_cells <= 0:
        raise ValueError("--min-detected-cells must be positive")

    n_select = min(int(n_variable_genes), MAX_VARIABLE_GENES, normalized.shape[0])
    means = np.asarray(normalized.mean(axis=1)).ravel()
    squared_means = np.asarray(normalized.power(2).mean(axis=1)).ravel()
    variances = squared_means - np.square(means)
    detected = np.asarray(normalized.getnnz(axis=1)).ravel()

    valid = np.flatnonzero(
        np.isfinite(variances)
        & (variances > 0)
        & (detected >= int(min_detected_cells))
    )
    if valid.size == 0:
        raise RuntimeError("No variable genes available after normalization and detection filtering")

    order = valid[np.lexsort((valid, -variances[valid]))]
    return order[: min(n_select, order.size)]


def split_group(indices: np.ndarray, max_cells: int, min_cells: int) -> list[np.ndarray]:
    if max_cells <= 0:
        raise ValueError("--max-cells-per-metacell must be positive")
    if min_cells <= 0:
        raise ValueError("--min-cells-per-metacell must be positive")

    if len(indices) <= max_cells:
        return [indices]

    n_chunks = int(np.ceil(len(indices) / float(max_cells)))
    while n_chunks > 1 and len(indices) // n_chunks < min_cells:
        n_chunks -= 1
    return [chunk.astype(int, copy=False) for chunk in np.array_split(indices, n_chunks) if len(chunk) > 0]


def subgroup_key(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    if not columns:
        return pd.Series(["all"] * len(frame), index=frame.index)
    return frame[columns].astype(str).agg("|".join, axis=1)


def build_metacell_groups(
    metadata: pd.DataFrame,
    min_cells: int,
    max_cells: int,
) -> tuple[list[np.ndarray], list[str], list[str]]:
    subgroup_columns = [column for column in SUBGROUP_COLUMNS if column in metadata.columns]
    groups: list[np.ndarray] = []
    labels: list[str] = []
    celltypes: list[str] = []

    for celltype, type_frame in metadata.groupby("celltype", sort=True):
        keys = subgroup_key(type_frame, subgroup_columns)
        base_groups: list[tuple[str, np.ndarray]] = []
        pooled_small: list[np.ndarray] = []

        for key, key_index in keys.groupby(keys, sort=True).groups.items():
            indices = np.asarray(type_frame.loc[key_index].index, dtype=int)
            if len(indices) >= min_cells:
                base_groups.append((str(key), indices))
            else:
                pooled_small.append(indices)

        if pooled_small:
            pooled = np.concatenate(pooled_small)
            if len(pooled) >= min_cells or not base_groups:
                base_groups.append(("pooled_small", pooled))

        if not base_groups:
            base_groups.append(("all", np.asarray(type_frame.index, dtype=int)))

        for key, indices in base_groups:
            for chunk_number, chunk in enumerate(split_group(np.sort(indices), max_cells, min_cells), start=1):
                groups.append(chunk)
                labels.append(f"{celltype}|{key}|chunk{chunk_number:02d}")
                celltypes.append(str(celltype))

    if len(groups) < 2 and len(metadata) >= 2:
        all_indices = np.arange(len(metadata), dtype=int)
        groups = [chunk.astype(int, copy=False) for chunk in np.array_split(all_indices, 2) if len(chunk) > 0]
        labels = [f"all|fallback|chunk{index + 1:02d}" for index in range(len(groups))]
        celltypes = ["all"] * len(groups)

    return groups, labels, celltypes


def build_metacell_profiles(
    expression: np.ndarray,
    metadata: pd.DataFrame,
    min_cells: int,
    max_cells: int,
) -> Metacells:
    groups, labels, celltypes = build_metacell_groups(metadata, min_cells, max_cells)
    if not groups:
        raise RuntimeError("No metacell groups could be constructed")

    profiles = np.vstack([expression[indices, :].mean(axis=0) for indices in groups]).astype(np.float32, copy=False)
    return Metacells(profiles=profiles, labels=labels, celltypes=celltypes)


def filter_metacell_variable_genes(
    expression: np.ndarray,
    gene_indices: np.ndarray,
    metacells: Metacells,
) -> tuple[np.ndarray, np.ndarray, Metacells]:
    if metacells.profiles.shape[0] < 2:
        return expression, gene_indices, metacells

    variances = metacells.profiles.var(axis=0)
    usable = np.flatnonzero(np.isfinite(variances) & (variances > 1e-8))
    if usable.size == 0:
        return expression, gene_indices, metacells

    filtered_metacells = Metacells(
        profiles=metacells.profiles[:, usable],
        labels=metacells.labels,
        celltypes=metacells.celltypes,
    )
    return expression[:, usable], gene_indices[usable], filtered_metacells


def compute_correlation_matrix(profiles: np.ndarray) -> np.ndarray:
    n_metacells, n_genes = profiles.shape
    if n_genes == 0:
        raise RuntimeError("No genes available for co-expression correlation")
    if n_genes == 1 or n_metacells < 2:
        return np.eye(n_genes, dtype=np.float32)

    centered = profiles.astype(np.float64, copy=False) - profiles.mean(axis=0, keepdims=True)
    norms = np.sqrt(np.sum(np.square(centered), axis=0))
    usable = norms > 1e-12
    standardized = np.zeros_like(centered)
    standardized[:, usable] = centered[:, usable] / norms[usable]

    correlation = standardized.T @ standardized
    correlation = np.nan_to_num(correlation, nan=0.0, posinf=0.0, neginf=0.0)
    correlation = (correlation + correlation.T) / 2.0
    correlation = np.clip(correlation, -1.0, 1.0)
    np.fill_diagonal(correlation, 1.0)
    return correlation.astype(np.float32, copy=False)


def requested_module_count(n_genes: int, min_module_size: int, max_modules: int) -> int:
    if min_module_size <= 0:
        raise ValueError("--min-module-size must be positive")
    if max_modules <= 0:
        raise ValueError("--max-modules must be positive")
    if n_genes <= min_module_size:
        return 1

    size_limited = max(1, n_genes // min_module_size)
    heuristic = max(2, int(round(np.sqrt(n_genes / 25.0))))
    return max(1, min(max_modules, size_limited, heuristic, n_genes))


def merge_small_modules(labels: np.ndarray, correlation: np.ndarray, min_module_size: int) -> np.ndarray:
    merged = labels.copy()
    while True:
        unique, counts = np.unique(merged, return_counts=True)
        sizes = dict(zip(unique.tolist(), counts.tolist()))
        small = [label for label, size in sizes.items() if size < min_module_size]
        large = [label for label, size in sizes.items() if size >= min_module_size]
        if not small or not large:
            return merged

        changed = False
        for label in sorted(small):
            gene_mask = merged == label
            best_label = None
            best_score = -np.inf
            for candidate in sorted(large):
                candidate_mask = merged == candidate
                score = float(np.nanmean(correlation[np.ix_(gene_mask, candidate_mask)]))
                if score > best_score:
                    best_score = score
                    best_label = candidate
            if best_label is not None:
                merged[gene_mask] = best_label
                changed = True

        if not changed:
            return merged


def cluster_gene_modules(
    correlation: np.ndarray,
    min_module_size: int,
    max_modules: int,
) -> np.ndarray:
    n_genes = correlation.shape[0]
    if n_genes == 0:
        raise RuntimeError("No genes available for module clustering")

    target_modules = requested_module_count(n_genes, min_module_size, max_modules)
    if target_modules == 1 or n_genes == 1:
        return np.ones(n_genes, dtype=int)

    distance = np.clip(1.0 - correlation.astype(np.float64, copy=False), 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)
    condensed = squareform(distance, checks=False)
    tree = linkage(condensed, method="average")
    labels = fcluster(tree, t=target_modules, criterion="maxclust").astype(int)
    return merge_small_modules(labels, correlation, min_module_size)


def name_modules(labels: np.ndarray) -> np.ndarray:
    unique = np.unique(labels)
    ordered = sorted(
        unique,
        key=lambda label: (-int(np.sum(labels == label)), int(np.flatnonzero(labels == label)[0])),
    )
    mapping = {label: f"module_{index:02d}" for index, label in enumerate(ordered, start=1)}
    return np.asarray([mapping[label] for label in labels], dtype=object)


def vector_correlations(matrix: np.ndarray, vector: np.ndarray) -> np.ndarray:
    if matrix.shape[0] < 2:
        return np.zeros(matrix.shape[1], dtype=np.float32)

    x = matrix.astype(np.float64, copy=False) - matrix.mean(axis=0, keepdims=True)
    y = vector.astype(np.float64, copy=False) - float(np.mean(vector))
    y_norm = np.sqrt(np.sum(np.square(y)))
    x_norm = np.sqrt(np.sum(np.square(x), axis=0))
    denom = x_norm * y_norm
    numerator = x.T @ y
    correlations = np.divide(numerator, denom, out=np.zeros_like(numerator), where=denom > 1e-12)
    return np.nan_to_num(correlations, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32, copy=False)


def first_eigengene(values: np.ndarray) -> np.ndarray:
    if values.shape[0] < 2:
        return np.zeros(values.shape[0], dtype=np.float32)

    centered = values.astype(np.float64, copy=False) - values.mean(axis=0, keepdims=True)
    scale = centered.std(axis=0)
    usable = scale > 1e-12
    if not np.any(usable):
        return np.zeros(values.shape[0], dtype=np.float32)

    z_values = np.zeros_like(centered)
    z_values[:, usable] = centered[:, usable] / scale[usable]
    u, singular_values, _ = np.linalg.svd(z_values, full_matrices=False)
    if singular_values.size == 0:
        return np.zeros(values.shape[0], dtype=np.float32)
    return (u[:, 0] * singular_values[0]).astype(np.float32, copy=False)


def compute_module_membership(
    profiles: np.ndarray,
    module_names: np.ndarray,
) -> np.ndarray:
    membership = np.zeros(profiles.shape[1], dtype=np.float32)
    for module in sorted(set(module_names.tolist())):
        indices = np.flatnonzero(module_names == module)
        module_values = profiles[:, indices]
        if len(indices) == 1:
            correlations = np.asarray([1.0], dtype=np.float32)
        else:
            eigengene = first_eigengene(module_values)
            correlations = vector_correlations(module_values, eigengene)
            mean_correlation = float(np.nanmean(correlations)) if correlations.size else 0.0
            if mean_correlation < 0:
                correlations = -correlations
        membership[indices] = correlations
    return membership


def zscore_columns(values: np.ndarray) -> np.ndarray:
    means = values.mean(axis=0, keepdims=True)
    sds = values.std(axis=0, keepdims=True)
    z_values = values.astype(np.float32, copy=True)
    z_values -= means.astype(np.float32, copy=False)
    np.divide(z_values, sds, out=z_values, where=sds > 1e-12)
    z_values[:, np.ravel(sds <= 1e-12)] = 0.0
    return z_values


def compute_module_scores(
    expression: np.ndarray,
    metadata: pd.DataFrame,
    module_names: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    modules = sorted(set(module_names.tolist()))
    z_expression = zscore_columns(expression)

    scores = pd.DataFrame({"barcode": metadata["barcode"].astype(str).to_numpy()})
    summary_source = pd.DataFrame({"celltype": metadata["celltype"].astype(str).to_numpy()})
    for module in modules:
        indices = np.flatnonzero(module_names == module)
        module_score = z_expression[:, indices].mean(axis=1) if len(indices) else np.zeros(len(metadata))
        scores[module] = module_score.astype(np.float32, copy=False)
        summary_source[module] = scores[module].to_numpy()

    counts = summary_source["celltype"].value_counts().rename("n_cells")
    means = summary_source.groupby("celltype", sort=True)[modules].mean()
    summary = means.join(counts).reset_index()
    summary = summary[["celltype", "n_cells", *modules]]
    return scores, summary, modules


def write_gene_modules(
    path: Path,
    genes: list[str],
    selected_indices: np.ndarray,
    expression: np.ndarray,
    module_names: np.ndarray,
    membership: np.ndarray,
) -> pd.DataFrame:
    module_sizes = pd.Series(module_names).value_counts().to_dict()
    result = pd.DataFrame(
        {
            "gene": [genes[index] for index in selected_indices],
            "gene_index": selected_indices.astype(int),
            "selected_rank": np.arange(1, len(selected_indices) + 1, dtype=int),
            "module": module_names,
            "module_size": [int(module_sizes[module]) for module in module_names],
            "module_membership_corr": membership,
            "mean_log_expression": expression.mean(axis=0),
            "variance_log_expression": expression.var(axis=0),
        }
    )
    result = result.sort_values(
        ["module", "module_membership_corr", "selected_rank"],
        ascending=[True, False, True],
    )
    result.to_csv(path, sep="\t", index=False)
    return result


def plot_module_celltype_heatmap(summary: pd.DataFrame, modules: list[str], path: Path) -> None:
    values = summary[modules].to_numpy(dtype=float) if modules else np.empty((len(summary), 0))
    if values.size == 0:
        values = np.zeros((max(1, len(summary)), 1), dtype=float)
        modules = ["module_score"]

    max_abs = float(np.nanmax(np.abs(values))) if np.isfinite(values).any() else 0.0
    limit = max(max_abs, 1e-6)

    width = max(5.0, 0.55 * len(modules) + 2.5)
    height = max(3.5, 0.38 * len(summary) + 2.0)
    fig, ax = plt.subplots(figsize=(width, height))
    image = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
    ax.set_xticks(np.arange(len(modules)))
    ax.set_xticklabels(modules, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(summary)))
    ax.set_yticklabels(summary["celltype"].astype(str).tolist())
    ax.set_xlabel("Co-expression module")
    ax.set_ylabel("Consensus cell type")
    ax.set_title("Mean module scores by cell type")

    if values.shape[0] * values.shape[1] <= 100:
        threshold = limit * 0.55
        for row in range(values.shape[0]):
            for column in range(values.shape[1]):
                value = values[row, column]
                color = "white" if abs(value) >= threshold else "black"
                ax.text(column, row, f"{value:.2f}", ha="center", va="center", color=color, fontsize=7)

    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Mean module score")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run(args: argparse.Namespace) -> str:
    np.random.seed(RANDOM_SEED)
    out_dir = args.out
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    paths = resolve_paths(args)
    require_files(paths)

    genes = read_vector(paths.genes, "genes.tsv")
    barcodes = read_barcodes(paths.barcodes)
    annotations, celltype_column = read_annotations(paths.annotations)
    counts = read_counts(paths.counts, genes, barcodes)
    counts, metadata = align_to_annotated_singlets(counts, annotations, barcodes)
    if counts.shape[1] != len(metadata):
        raise RuntimeError("Aligned count matrix and metadata row counts differ")
    if counts.shape[1] < 2:
        raise RuntimeError("At least two annotated singlets are required for co-expression analysis")

    normalized = normalize_log1p_gene_by_cell(counts, args.normalization_scale)
    selected_indices = select_variable_genes(normalized, args.n_variable_genes, args.min_detected_cells)
    expression = normalized[selected_indices, :].T.toarray().astype(np.float32, copy=False)

    metacells = build_metacell_profiles(
        expression=expression,
        metadata=metadata,
        min_cells=args.min_cells_per_metacell,
        max_cells=args.max_cells_per_metacell,
    )
    expression, selected_indices, metacells = filter_metacell_variable_genes(expression, selected_indices, metacells)
    if len(selected_indices) == 0:
        raise RuntimeError("No selected genes retained enough metacell-level variation")

    correlation = compute_correlation_matrix(metacells.profiles)
    raw_labels = cluster_gene_modules(correlation, args.min_module_size, args.max_modules)
    module_names = name_modules(raw_labels)
    membership = compute_module_membership(metacells.profiles, module_names)

    gene_modules_path = tables_dir / "coexpression_gene_modules.tsv"
    module_scores_path = tables_dir / "coexpression_module_scores.tsv"
    summary_path = tables_dir / "coexpression_module_celltype_summary.tsv"
    heatmap_path = figures_dir / "coexpression_module_celltype_heatmap.png"

    write_gene_modules(gene_modules_path, genes, selected_indices, expression, module_names, membership)
    scores, summary, modules = compute_module_scores(expression, metadata, module_names)
    scores.to_csv(module_scores_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    plot_module_celltype_heatmap(summary, modules, heatmap_path)

    detail = (
        f"cells={counts.shape[1]} selected_genes={len(selected_indices)} modules={len(modules)} "
        f"metacells={metacells.profiles.shape[0]} celltypes={metadata['celltype'].nunique()} "
        f"celltype_column={celltype_column}"
    )
    output_summary = (
        f"{detail}; outputs={gene_modules_path}, {module_scores_path}, "
        f"{summary_path}, {heatmap_path}"
    )
    append_status(out_dir, MODULE, "PASS", detail)
    return output_summary


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out_dir = args.out
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    append_status(out_dir, MODULE, "START", "local_hdWGCNA_style_coexpression")
    try:
        summary = run(args)
    except Exception as exc:
        try:
            append_status(out_dir, MODULE, "FAIL", f"{type(exc).__name__}: {exc}")
        except Exception:
            pass
        raise
    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
