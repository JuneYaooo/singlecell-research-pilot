#!/usr/bin/env python3
"""Validate Seurat V5 course-code integration for this skill."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from zipfile import ZipFile


REQUIRED_REFS = [
    Path("references/singlecell-full-workflow-user-facing.md"),
    Path("references/seurat-v5-workflows.md"),
    Path("references/seurat-v5-code-index.md"),
]

REQUIRED_PUBLIC_FILES = [
    Path("README.md"),
    Path(".gitignore"),
]

PUBLIC_NO_LOCAL_PATH_FILES = [
    Path("SKILL.md"),
    Path("README.md"),
    Path("references/workflows.md"),
    Path("references/singlecell-full-workflow-user-facing.md"),
    Path("references/seurat-v5-workflows.md"),
    Path("references/seurat-v5-code-index.md"),
    Path("scripts/run_full_course_logic.py"),
    Path("scripts/run_logic_phase2_all.sh"),
    Path("scripts/run_seurat_v5_course_core.R"),
    Path("scripts/run_seurat_v5_course_cibersort.R"),
]

FORBIDDEN_PUBLIC_PATH_FRAGMENTS = [
    "/" + "Users" + "/june2",
    "/" + "Users" + "/",
    "code" + "/data",
    "github" + "/singlecell-research-pilot",
]

REQUIRED_WORKFLOW_TERMS = [
    "Seurat V5 Course-Derived Workflow Reference",
    "Input Import",
    "QC and Filtering",
    "Normalization, Integration, and Clustering",
    "Annotation",
    "Downstream Modules",
    "Execution Guardrails",
]

REQUIRED_USER_FACING_TERMS = [
    "Single-cell Full Workflow User-facing Reference",
    "Data Import",
    "Quality Control",
    "Clustering",
    "Cell Annotation",
    "Cell-cell Communication",
    "Pseudotime",
    "CytoTRACE",
    "Subset Extraction",
    "FindMarkers",
    "FindAllMarkers",
    "Gene Expression Plotting",
    "File Management",
    "Validation Matrix",
    "Replacement Policy",
]


def decode_zip_name(name: str) -> str:
    """Decode legacy Windows ZIP names that were stored without UTF-8 flags."""
    try:
        return name.encode("cp437").decode("gbk")
    except UnicodeError:
        return name


def list_r_scripts(zip_path: Path) -> list[str]:
    with ZipFile(zip_path) as archive:
        scripts = []
        for info in archive.infolist():
            decoded = decode_zip_name(info.filename)
            if decoded.lower().endswith(".r"):
                scripts.append(decoded)
        return sorted(scripts)


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def validate(skill_root: Path, zip_path: Path | None) -> list[str]:
    errors: list[str] = []

    skill_md = skill_root / "SKILL.md"
    require(skill_md.exists(), "SKILL.md is missing", errors)
    skill_text = skill_md.read_text(encoding="utf-8") if skill_md.exists() else ""

    for public_file in REQUIRED_PUBLIC_FILES:
        require((skill_root / public_file).exists(), f"{public_file} is missing", errors)

    for public_file in PUBLIC_NO_LOCAL_PATH_FILES:
        path = skill_root / public_file
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        for fragment in FORBIDDEN_PUBLIC_PATH_FRAGMENTS:
            require(fragment not in text, f"{public_file} contains local path fragment: {fragment}", errors)

    for ref in REQUIRED_REFS:
        require((skill_root / ref).exists(), f"{ref} is missing", errors)
        require(str(ref) in skill_text, f"SKILL.md does not link to {ref}", errors)

    workflow_path = skill_root / "references/seurat-v5-workflows.md"
    workflow_text = workflow_path.read_text(encoding="utf-8") if workflow_path.exists() else ""
    for term in REQUIRED_WORKFLOW_TERMS:
        require(term in workflow_text, f"workflow reference missing required term: {term}", errors)

    user_facing_path = skill_root / "references/singlecell-full-workflow-user-facing.md"
    user_facing_text = user_facing_path.read_text(encoding="utf-8") if user_facing_path.exists() else ""
    for term in REQUIRED_USER_FACING_TERMS:
        require(term in user_facing_text, f"user-facing workflow reference missing required term: {term}", errors)

    index_path = skill_root / "references/seurat-v5-code-index.md"
    index_text = index_path.read_text(encoding="utf-8") if index_path.exists() else ""
    require("Total R scripts: 35" in index_text, "code index does not record 35 R scripts", errors)

    if zip_path is not None:
        require(zip_path.exists(), f"zip path does not exist: {zip_path}", errors)
        if zip_path.exists():
            scripts = list_r_scripts(zip_path)
            require(len(scripts) == 35, f"expected 35 R scripts in zip, found {len(scripts)}", errors)
            for script in scripts:
                require(script in index_text, f"code index missing script path from zip: {script}", errors)

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Skill root directory. Defaults to the parent of scripts/.",
    )
    parser.add_argument(
        "--zip",
        type=Path,
        default=None,
        help="Optional Seurat V5 course ZIP to validate script coverage against.",
    )
    args = parser.parse_args()

    errors = validate(args.root.resolve(), args.zip.resolve() if args.zip else None)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print("Seurat V5 skill integration is valid.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
