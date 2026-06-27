#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-analysis/seurat_v5_logic_run}"
PY="${ROOT_DIR}/.conda/scverse-course/bin/python"
RSCRIPT="${ROOT_DIR}/.conda/seurat-core/bin/Rscript"
COURSE_ROOT="${COURSE_ROOT:-${SEURAT_V5_COURSE_ROOT:-}}"

cd "${ROOT_DIR}"

if [[ ! -x "${PY}" ]]; then
  echo "Python environment not found: ${PY}" >&2
  exit 1
fi
if [[ ! -x "${RSCRIPT}" ]]; then
  echo "Rscript environment not found: ${RSCRIPT}" >&2
  exit 1
fi
if [[ -z "${COURSE_ROOT}" ]]; then
  echo "COURSE_ROOT or SEURAT_V5_COURSE_ROOT is required for source traceability" >&2
  exit 1
fi

if [[ ! -s "${OUT_DIR}/tables/consensus_annotation.tsv" ]]; then
  echo "Phase 1 consensus_annotation.tsv missing under ${OUT_DIR}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}/logs"

run_module() {
  local name="$1"
  shift
  echo "=== Running ${name} ==="
  "$@" > "${OUT_DIR}/logs/${name}.log" 2>&1
  tail -n 20 "${OUT_DIR}/logs/${name}.log"
}

run_module export_consensus_embeddings "${RSCRIPT}" scripts/run_logic_export_consensus_embeddings.R --out "${OUT_DIR}"
run_module communication "${PY}" scripts/run_logic_communication.py --out "${OUT_DIR}"
run_module trajectory "${PY}" scripts/run_logic_trajectory.py --out "${OUT_DIR}"
run_module cytotrace "${PY}" scripts/run_logic_cytotrace.py --out "${OUT_DIR}"
run_module subset_extraction "${RSCRIPT}" scripts/run_logic_subset.R --out "${OUT_DIR}"
run_module findmarkers_pairwise "${RSCRIPT}" scripts/run_logic_findmarkers.R --out "${OUT_DIR}"
run_module gene_expression_plotting "${RSCRIPT}" scripts/run_logic_gene_expression_plotting.R --out "${OUT_DIR}"
run_module cnv_proxy "${PY}" scripts/run_logic_cnv_proxy.py --out "${OUT_DIR}"
run_module coexpression "${PY}" scripts/run_logic_coexpression.py --out "${OUT_DIR}"
run_module celltypist_bounded "${PY}" scripts/run_logic_celltypist_bounded.py --out "${OUT_DIR}"
run_module source_traceability "${PY}" scripts/build_source_traceability.py --course-root "${COURSE_ROOT}" --logic-dir "${OUT_DIR}"
run_module file_manifest "${PY}" scripts/run_logic_file_manifest.py --out "${OUT_DIR}"
run_module final_report "${PY}" scripts/run_logic_final_report.py --out "${OUT_DIR}"
run_module final_artifact_manifest "${PY}" scripts/run_logic_file_manifest.py --out "${OUT_DIR}" --final

echo "Phase 2 run complete. See ${OUT_DIR}/final_report.md"
