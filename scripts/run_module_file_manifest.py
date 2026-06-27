#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


MODULE = "file_manifest"
FINAL_MODULE = "final_artifact_manifest"
EXCLUDED_NAMES = {"file_manifest.tsv", "file_manifest_summary.tsv", "final_artifact_manifest.tsv", "final_artifact_manifest_summary.tsv"}


def append_status(out_dir: Path, module: str, status: str, detail: str = "") -> None:
    status_file = out_dir / "tables" / "module_status.tsv"
    status_file.parent.mkdir(parents=True, exist_ok=True)
    if not status_file.exists():
        status_file.write_text("module\tstatus\tdetail\n")
    clean_detail = " ".join(str(detail).replace("\t", " ").splitlines())
    with status_file.open("a") as handle:
        handle.write(f"{module}\t{status}\t{clean_detail}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a reproducibility/download manifest for a logic run directory.")
    parser.add_argument("--out", required=True, type=Path, help="Logic output directory.")
    parser.add_argument("--final", action="store_true", help="Write final_artifact_manifest.tsv and include manifest files.")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def module_from_path(rel_path: str) -> str:
    parts = rel_path.split("/")
    name = parts[-1]
    if rel_path.startswith("figures/cytotrace") or rel_path.startswith("tables/cytotrace"):
        return "cytotrace"
    if name.startswith("annotated_consensus_"):
        return "export_consensus_embeddings"
    if rel_path.startswith("objects/rds/subset_") or name.startswith("subset_"):
        return "subset_extraction"
    if name.startswith("findmarkers_pairwise"):
        return "findmarkers_pairwise"
    if name.startswith("gene_expression_"):
        return "gene_expression_plotting"
    if name.startswith("communication_"):
        return "communication_lr"
    if name.startswith("trajectory_"):
        return "trajectory"
    if name.startswith("cnv_proxy"):
        return "cnv_proxy"
    if name.startswith("coexpression_"):
        return "coexpression"
    if name.startswith("celltypist_"):
        return "celltypist_bounded"
    if "source_traceability" in name:
        return "source_traceability"
    if name.startswith("final_") or name in {"final_report.md"}:
        if name.startswith("final_artifact_manifest"):
            return "final_artifact_manifest"
        return "final_report"
    if name.startswith("module_status"):
        return "status"
    if rel_path.startswith("logs/"):
        return "logs"
    if rel_path.startswith("objects/export/"):
        return "export_seurat"
    return "general"


def file_kind(path: Path) -> str:
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".tsv") or suffixes.endswith(".csv"):
        return "table"
    if suffixes.endswith(".png") or suffixes.endswith(".pdf") or suffixes.endswith(".svg"):
        return "figure"
    if suffixes.endswith(".rds") or suffixes.endswith(".h5ad") or suffixes.endswith(".mtx.gz"):
        return "object"
    if suffixes.endswith(".log") or suffixes.endswith(".txt"):
        return "log"
    if suffixes.endswith(".md") or suffixes.endswith(".html"):
        return "report"
    return "other"


def build_manifest(out_dir: Path, include_self: bool = False) -> pd.DataFrame:
    rows = []
    for path in sorted(out_dir.rglob("*")):
        if not path.is_file():
            continue
        if not include_self and path.name in EXCLUDED_NAMES:
            continue
        rel_path = path.relative_to(out_dir).as_posix()
        rows.append(
            {
                "path": rel_path,
                "module": module_from_path(rel_path),
                "kind": file_kind(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    if not rows:
        raise RuntimeError(f"No files found under {out_dir}")
    return pd.DataFrame(rows).sort_values(["kind", "module", "path"]).reset_index(drop=True)


def latest_status(status: pd.DataFrame) -> pd.DataFrame:
    if status.empty:
        return pd.DataFrame(columns=["module", "status", "detail"])
    rows = []
    for module, group in status.groupby("module", sort=False):
        non_start = group[group["status"] != "START"]
        rows.append((non_start if not non_start.empty else group).iloc[-1])
    return pd.DataFrame(rows).reset_index(drop=True)


def write_final_module_status(out_dir: Path) -> None:
    status_path = out_dir / "tables" / "module_status.tsv"
    if not status_path.exists():
        raise RuntimeError(f"Missing module status table: {status_path}")
    status_events = pd.read_csv(status_path, sep="\t")
    latest_status(status_events).to_csv(out_dir / "tables" / "final_module_status.tsv", sep="\t", index=False)


def run(out_dir: Path, *, final: bool = False) -> str:
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)
    module = FINAL_MODULE if final else MODULE
    output_name = "final_artifact_manifest.tsv" if final else "file_manifest.tsv"
    summary_name = "final_artifact_manifest_summary.tsv" if final else "file_manifest_summary.tsv"
    append_status(out_dir, module, "START", f"building {output_name}")

    manifest_path = out_dir / "tables" / output_name
    summary_path = out_dir / "tables" / summary_name
    if final and not manifest_path.exists():
        manifest_path.write_text("path\tmodule\tkind\tsize_bytes\tsha256\n")
    if final and not summary_path.exists():
        summary_path.write_text("kind\tmodule\tfiles\ttotal_bytes\n")

    manifest = build_manifest(out_dir, include_self=final)
    manifest.to_csv(manifest_path, sep="\t", index=False)

    summary = (
        manifest.groupby(["kind", "module"], dropna=False)
        .agg(files=("path", "count"), total_bytes=("size_bytes", "sum"))
        .reset_index()
        .sort_values(["kind", "module"])
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    detail = f"files={len(manifest)} total_bytes={int(manifest['size_bytes'].sum())} sha256=yes"
    append_status(out_dir, module, "PASS", detail)
    if final:
        write_final_module_status(out_dir)
    return detail


def main() -> None:
    args = parse_args()
    try:
        detail = run(args.out, final=args.final)
    except Exception as exc:
        try:
            append_status(args.out, FINAL_MODULE if args.final else MODULE, "FAIL", f"{type(exc).__name__}:{exc}")
        except Exception:
            pass
        raise
    print(detail)


if __name__ == "__main__":
    main()
