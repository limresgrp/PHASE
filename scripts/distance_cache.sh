#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/offline_select.sh"
OFFLINE_ROOT=""
OFFLINE_PROJECT_ID=""
OFFLINE_SYSTEM_ID=""
CLUSTER_ID=""
while [ "$#" -gt 0 ]; do
  case "$1" in
    --root) OFFLINE_ROOT="$2"; shift 2 ;;
    --project-id) OFFLINE_PROJECT_ID="$2"; shift 2 ;;
    --system-id) OFFLINE_SYSTEM_ID="$2"; shift 2 ;;
    --cluster-id) CLUSTER_ID="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
if [ -z "${VIRTUAL_ENV:-}" ] || [ ! -x "${VIRTUAL_ENV}/bin/python" ]; then
  echo "Activate the PHASE virtual environment first." >&2
  exit 1
fi
if [ -z "$OFFLINE_ROOT" ]; then offline_prompt_root "${PHASE_DATA_ROOT:-${ROOT_DIR}/data}"; fi
export PHASE_DATA_ROOT="$OFFLINE_ROOT"
if [ -z "$OFFLINE_PROJECT_ID" ]; then offline_select_project; fi
if [ -z "$OFFLINE_SYSTEM_ID" ]; then offline_select_system; fi
if [ -z "$CLUSTER_ID" ]; then
  CLUSTER_ROW="$(offline_select_cluster)"
  CLUSTER_ID="${CLUSTER_ROW%%|*}"
fi
[ -n "$CLUSTER_ID" ] || exit 1
BASE=("${VIRTUAL_ENV}/bin/python" -m phase.scripts.distance_cache --root "$OFFLINE_ROOT" --project-id "$OFFLINE_PROJECT_ID" --system-id "$OFFLINE_SYSTEM_ID" --cluster-id "$CLUSTER_ID")
ACTION_ROW="$(offline_choose_one "Distance analysis:" $'build|Build/refresh aligned CA distances\nstatus|Validate cache / inspect metadata\npair|Inspect pair distance histogram\njs|Compute conditional pair JS')"
case "${ACTION_ROW%%|*}" in
  build)
    CUTOFF="$(prompt "Candidate-edge cutoff in angstrom (independent of Potts cutoff)" "12")"
    CMD=("${BASE[@]}" build --max-distance-cutoff "$CUTOFF" --progress)
    if ! prompt_bool "Use minimum-image distances when a valid periodic box exists? (Y/n)" "Y"; then CMD+=(--no-pbc); fi
    if prompt_bool "Override stored trajectory paths? (y/N)" "N"; then
      while true; do
        OVERRIDE="$(prompt "STATE_ID=/absolute/trajectory/path (blank to finish)" "")"
        [ -n "$OVERRIDE" ] || break
        CMD+=(--trajectory "$OVERRIDE")
      done
    fi
    ;;
  js)
    STATE_LINES="$(_offline_list list-states --project-id "$OFFLINE_PROJECT_ID" --system-id "$OFFLINE_SYSTEM_ID")"
    A="$(offline_choose_multi "Condition A states (must have merged cluster frames):" "$STATE_LINES")"
    B="$(offline_choose_multi "Condition B states (must have merged cluster frames):" "$STATE_LINES")"
    A_IDS="$(printf '%s\n' "$A" | awk -F'|' 'NF {print $1}' | paste -sd, -)"
    B_IDS="$(printf '%s\n' "$B" | awk -F'|' 'NF {print $1}' | paste -sd, -)"
    [ -n "$A_IDS" ] && [ -n "$B_IDS" ] || exit 1
    THRESHOLDS="$(prompt "Distance thresholds, comma separated (blank=global)" "")"
    MINIMUM="$(prompt "Minimum valid frames in each condition per bin" "1")"
    NAME="$(prompt "Analysis name (blank=automatic)" "")"
    META_A="$(prompt "Optional metastable IDs in condition A, comma separated" "")"
    META_B="$(prompt "Optional metastable IDs in condition B, comma separated" "")"
    CMD=("${BASE[@]}" js --states-a "$A_IDS" --states-b "$B_IDS" --thresholds "$THRESHOLDS" --min-bin-frames "$MINIMUM")
    if [ -n "$NAME" ]; then CMD+=(--name "$NAME"); fi
    if [ -n "$META_A" ]; then CMD+=(--metastable-a "$META_A"); fi
    if [ -n "$META_B" ]; then CMD+=(--metastable-b "$META_B"); fi
    ;;
  status) CMD=("${BASE[@]}" status) ;;
  pair)
    R="$(prompt "First residue key or zero-based column index" "res_2")"
    S="$(prompt "Second residue key or zero-based column index" "res_3")"
    STATES="$(prompt "Macro state IDs (comma separated, blank=all)" "")"
    METAS="$(prompt "Metastable IDs (comma separated, blank=all)" "")"
    THRESHOLDS="$(prompt "Distance thresholds (comma separated, blank=one bin)" "")"
    CMD=("${BASE[@]}" pair --pair "$R" "$S" --thresholds "$THRESHOLDS")
    if [ -n "$STATES" ]; then CMD+=(--state-ids "$STATES"); fi
    if [ -n "$METAS" ]; then CMD+=(--metastable-ids "$METAS"); fi
    ;;
  *) exit 0 ;;
esac
"${CMD[@]}"
