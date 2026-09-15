#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/phase_defaults.sh"
source "${ROOT_DIR}/scripts/offline_select.sh"

if [ -n "${VIRTUAL_ENV:-}" ] && [ -x "${VIRTUAL_ENV}/bin/python" ]; then
  PYTHON_BIN="${VIRTUAL_ENV}/bin/python"
else
  echo "Activate a PHASE virtual environment first (for example: source .venv-phase/bin/activate)." >&2
  exit 1
fi

OFFLINE_ROOT="${PHASE_DATA_ROOT:-${ROOT_DIR}/data}"
offline_prompt_root "$OFFLINE_ROOT"
offline_select_project
offline_select_system
CLUSTER_ROW="$(offline_select_cluster)"
CLUSTER_ID="$(printf "%s" "$CLUSTER_ROW" | awk -F'|' '{print $1}')"
CLUSTER_PATH="$(printf "%s" "$CLUSTER_ROW" | awk -F'|' '{print $3}')"
[ -z "$CLUSTER_ID" ] && echo "No cluster selected." >&2 && exit 1

STATE_LINES="$(_offline_list list-states --project-id "$OFFLINE_PROJECT_ID" --system-id "$OFFLINE_SYSTEM_ID")"
[ -z "$(trim "$STATE_LINES")" ] && echo "No states found." >&2 && exit 1

read_state_selection() {
  local topology_path="$1"
  local metadata_path
  metadata_path="$(dirname "$topology_path")/state_metadata.json"
  if [ ! -f "$metadata_path" ]; then
    printf "protein"
    return
  fi
  "$PYTHON_BIN" -c 'import json, sys; print(str(json.load(open(sys.argv[1])).get("residue_selection") or "protein"))' "$metadata_path"
}

read_state_resid_shift() {
  local topology_path="$1"
  local metadata_path
  metadata_path="$(dirname "$topology_path")/state_metadata.json"
  if [ ! -f "$metadata_path" ]; then
    printf "0"
    return
  fi
  "$PYTHON_BIN" -c 'import json, sys; print(int(json.load(open(sys.argv[1])).get("resid_shift") or 0))' "$metadata_path"
}

state_spec() {
  local row="$1"
  local state_id name pdb traj selection resid_shift
  state_id="$(printf "%s" "$row" | awk -F'|' '{print $1}')"
  name="$(printf "%s" "$row" | awk -F'|' '{print $2}')"
  pdb="$(printf "%s" "$row" | awk -F'|' '{print $3}')"
  traj="$(printf "%s" "$row" | awk -F'|' '{print $4}')"
  local default_path="${traj:-$pdb}"
  local system_root="${OFFLINE_ROOT}/projects/${OFFLINE_PROJECT_ID}/systems/${OFFLINE_SYSTEM_ID}"
  [[ "$default_path" = /* ]] || default_path="${system_root}/${default_path}"
  [[ "$pdb" = /* ]] || pdb="${system_root}/${pdb}"
  local selected_path
  selected_path="$(prompt "Trajectory for ${name} (blank uses ${default_path})" "$default_path")"
  selected_path="$(trim "$selected_path")"
  [ -z "$selected_path" ] && selected_path="$default_path"
  [[ "$selected_path" = /* ]] || selected_path="$(pwd)/${selected_path}"
  [ -z "$selected_path" ] && echo "No trajectory supplied for ${name}." >&2 && return 1
  selection="$(read_state_selection "$pdb")"
  resid_shift="$(read_state_resid_shift "$pdb")"
  printf "%s|%s|%s|%s|%s" "$name" "$selected_path" "$pdb" "$selection" "$resid_shift"
}

ENDPOINT_A_ROW="$(offline_choose_one "Select endpoint A state:" "$STATE_LINES")"
ENDPOINT_B_ROW="$(offline_choose_one "Select endpoint B state:" "$STATE_LINES")"
ENDPOINT_A="$(state_spec "$ENDPOINT_A_ROW")"
ENDPOINT_B="$(state_spec "$ENDPOINT_B_ROW")"
[ "$(printf "%s" "$ENDPOINT_A_ROW" | awk -F'|' '{print $1}')" = "$(printf "%s" "$ENDPOINT_B_ROW" | awk -F'|' '{print $1}')" ] && echo "Endpoints must be different states." >&2 && exit 1

TRAJ_LINES="$(_offline_list list-trajectories --project-id "$OFFLINE_PROJECT_ID" --system-id "$OFFLINE_SYSTEM_ID")"
[ -z "$(trim "$TRAJ_LINES")" ] && echo "No state trajectories found." >&2 && exit 1
TRAJ_ROWS="$(offline_choose_multi "Select trajectories to evaluate:" "$TRAJ_LINES")"
[ -z "$(trim "$TRAJ_ROWS")" ] && echo "No trajectories selected." >&2 && exit 1

TRAJ_ARGS=()
# Keep the selected rows on descriptor 3 so each trajectory override prompt
# still reads from the user's terminal on stdin.
while IFS= read -r -u 3 row; do
  [ -z "$(trim "$row")" ] && continue
  state_id="$(printf "%s" "$row" | awk -F'|' '{print $1}')"
  state_name="$(printf "%s" "$row" | awk -F'|' '{print $2}')"
  trajectory_path="$(printf "%s" "$row" | awk -F'|' '{print $3}')"
  state_row="$(printf "%s\n" "$STATE_LINES" | awk -F'|' -v sid="$state_id" '$1==sid {print; exit}')"
  topology_path="$(printf "%s" "$state_row" | awk -F'|' '{print $3}')"
  [ -z "$state_row" ] && echo "State not found for trajectory: $state_id" >&2 && exit 1
  system_root="${OFFLINE_ROOT}/projects/${OFFLINE_PROJECT_ID}/systems/${OFFLINE_SYSTEM_ID}"
  [[ "$trajectory_path" = /* ]] || trajectory_path="${system_root}/${trajectory_path}"
  [[ "$topology_path" = /* ]] || topology_path="${system_root}/${topology_path}"
  trajectory_spec="$(prompt "Trajectory for ${state_name} (blank uses ${trajectory_path})" "$trajectory_path")"
  trajectory_spec="$(trim "$trajectory_spec")"
  [ -z "$trajectory_spec" ] && trajectory_spec="$trajectory_path"
  [[ "$trajectory_spec" = /* ]] || trajectory_spec="$(pwd)/${trajectory_spec}"
  selection="$(read_state_selection "$topology_path")"
  resid_shift="$(read_state_resid_shift "$topology_path")"
  trajectory_spec="${state_name}|${trajectory_spec}|${topology_path}|${selection}|${resid_shift}"
  TRAJ_ARGS+=(--trajectory "$trajectory_spec")
done 3<<< "$TRAJ_ROWS"

CLUSTER_DIR="$(dirname "${CLUSTER_PATH:-${OFFLINE_ROOT}/projects/${OFFLINE_PROJECT_ID}/systems/${OFFLINE_SYSTEM_ID}/clusters/${CLUSTER_ID}/cluster.npz}")"
DEFAULT_OUTPUT="${CLUSTER_DIR}/analyses/rmsd_ordering"
OUTPUT_DIR="$(prompt "Output directory" "$DEFAULT_OUTPUT")"
ALIGN="$(prompt "Align structures before RMSD? (Y/n)" "Y")"
STRIDE="$(prompt "Frame stride for evaluated trajectories" "100")"
WORKERS="$(prompt "Parallel trajectory workers" "1")"

ENDPOINT_A_SPEC="$ENDPOINT_A"
ENDPOINT_B_SPEC="$ENDPOINT_B"
CMD=("$PYTHON_BIN" "${ROOT_DIR}/scripts/rmsd_ordering.py" --endpoint-a "$ENDPOINT_A_SPEC" --endpoint-b "$ENDPOINT_B_SPEC" --output-dir "$OUTPUT_DIR" --stride "$STRIDE" --workers "$WORKERS")
CMD+=("${TRAJ_ARGS[@]}")
[ "$(printf "%s" "$ALIGN" | tr '[:upper:]' '[:lower:]')" = "n" ] && CMD+=(--no-alignment)

echo ""
printf 'Running: '; printf '%q ' "${CMD[@]}"; echo
"${CMD[@]}"
