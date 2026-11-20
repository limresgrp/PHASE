#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/offline_select.sh"

DEFAULT_ROOT="${PHASE_DATA_ROOT:-${ROOT_DIR}/data}"

prompt() {
  local label="$1"
  local default="$2"
  local var
  read -r -p "${label} [${default}]: " var
  if [ -z "$var" ]; then
    echo "$default"
  else
    echo "$var"
  fi
}

trim() {
  printf "%s" "$1" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
}

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

if [ -n "${VIRTUAL_ENV:-}" ] && [ -x "${VIRTUAL_ENV}/bin/python" ]; then
  PYTHON_BIN="${VIRTUAL_ENV}/bin/python"
else
  echo "No active virtual environment detected." >&2
  echo "Activate .venv-phase or .venv-potts-fit first." >&2
  exit 1
fi

if [ -z "$OFFLINE_ROOT" ]; then
  offline_prompt_root "$DEFAULT_ROOT"
else
  OFFLINE_ROOT="$(trim "$OFFLINE_ROOT")"
  export PHASE_DATA_ROOT="$OFFLINE_ROOT"
fi
if [ -z "$OFFLINE_PROJECT_ID" ]; then offline_select_project; fi
if [ -z "$OFFLINE_SYSTEM_ID" ]; then offline_select_system; fi
if [ -z "$CLUSTER_ID" ]; then
  CLUSTER_ROW="$(offline_select_cluster)"
  CLUSTER_ID="$(printf "%s" "$CLUSTER_ROW" | awk -F'|' '{print $1}')"
fi
if [ -z "$CLUSTER_ID" ]; then
  echo "No cluster selected." >&2
  exit 1
fi

MODEL_LINES="$(_offline_list list-models --project-id "$OFFLINE_PROJECT_ID" --system-id "$OFFLINE_SYSTEM_ID" || true)"
MODEL_LINES="$(printf "%s\n" "$MODEL_LINES" | awk -F'|' -v cid="$CLUSTER_ID" '$4==cid')"
if [ -z "$(trim "$MODEL_LINES")" ]; then
  echo "No Potts models found for cluster: $CLUSTER_ID" >&2
  exit 1
fi

MODEL_ARGS=()
BATCH_MODE=false
if prompt_bool "Advanced: run all unique pairs of multiple models? (y/N)" "N"; then
  BATCH_MODE=true
  MODEL_ROWS="$(offline_choose_multi "Select models (each unique pair creates an independent analysis):" "$MODEL_LINES")"
  MODEL_IDS="$(printf "%s\n" "$MODEL_ROWS" | awk -F'|' 'NF {print $1}' | awk '!seen[$0]++' | paste -sd',' -)"
  [[ "$MODEL_IDS" == *,* ]] || { echo "Select at least two models." >&2; exit 1; }
  MODEL_ARGS=(--model-ids "$MODEL_IDS")
  echo "Pair orientation follows this model order: $MODEL_IDS"
else
  MODEL_A_ROW="$(offline_choose_one "Select model A:" "$MODEL_LINES")"
  MODEL_A_ID="$(printf "%s" "$MODEL_A_ROW" | awk -F'|' '{print $1}')"
  MODEL_B_ROW="$(offline_choose_one "Select model B:" "$MODEL_LINES")"
  MODEL_B_ID="$(printf "%s" "$MODEL_B_ROW" | awk -F'|' '{print $1}')"
  if [ -z "$MODEL_A_ID" ] || [ -z "$MODEL_B_ID" ]; then
    echo "Model A and model B are required." >&2
    exit 1
  fi
  if [ "$MODEL_A_ID" = "$MODEL_B_ID" ]; then
    echo "Model A and model B must be different." >&2
    exit 1
  fi
  MODEL_ARGS=(--model-a-id "$MODEL_A_ID" --model-b-id "$MODEL_B_ID")
fi

MODE="$(prompt "Mode (ensemble/parameter)" "ensemble")"
CMD=("$PYTHON_BIN" -m phase.scripts.potts_model_pair_hotspots
  --root "$OFFLINE_ROOT" --project-id "$OFFLINE_PROJECT_ID" --system-id "$OFFLINE_SYSTEM_ID"
  --cluster-id "$CLUSTER_ID" "${MODEL_ARGS[@]}" --mode "$MODE")
if [ "$MODE" = "ensemble" ]; then
  echo "Select assigned MD reference samples. Use independent endpoint trajectories; held-out replicas are preferable."
  SAMPLE_LINES="$(_offline_list list-sampling --project-id "$OFFLINE_PROJECT_ID" --system-id "$OFFLINE_SYSTEM_ID")"
  ALL_SAMPLE_LINES="$(printf "%s\n" "$SAMPLE_LINES" | awk -F'|' -v cid="$CLUSTER_ID" '$1==cid')"
  SAMPLE_LINES="$(printf "%s\n" "$SAMPLE_LINES" | awk -F'|' -v cid="$CLUSTER_ID" '$1==cid && $3=="md_eval"')"
  # The Python validator rejects non-MD reference choices, overlapping frames and raw delta models.
  if [ "$BATCH_MODE" = true ]; then
    while IFS='|' read -r -u 3 model_id model_name _; do
      [ -z "$model_id" ] && continue
      ROWS="$(offline_choose_multi "Select endpoint MD references for ${model_name} (${model_id}):" "$SAMPLE_LINES")"
      if [ -z "$(trim "$ROWS")" ]; then echo "Select references for every model." >&2; exit 1; fi
      IDS="$(printf "%s\n" "$ROWS" | awk -F'|' '{print $4}' | while IFS= read -r p; do basename "$(dirname "$p")"; done | paste -sd',' -)"
      CMD+=(--model-reference "${model_id}:${IDS}")
    done 3<<< "$MODEL_ROWS"
  else
    for GROUP in a b; do
      ROWS="$(offline_choose_multi "Select endpoint ${GROUP} MD reference sample(s):" "$SAMPLE_LINES")"
      if [ -z "$(trim "$ROWS")" ]; then echo "Select at least one reference MD sample for each endpoint." >&2; exit 1; fi
      IDS="$(printf "%s\n" "$ROWS" | awk -F'|' '{print $4}' | while IFS= read -r p; do basename "$(dirname "$p")"; done | paste -sd',' -)"
      CMD+=("--reference-${GROUP}" "$IDS")
    done
  fi
  if prompt_bool "Include additional ensembles for preference coloring? (y/N)" "N"; then
    ROWS="$(offline_choose_multi "Select additional MD or sampled ensembles:" "$ALL_SAMPLE_LINES")"
    if [ -n "$(trim "$ROWS")" ]; then
      IDS="$(printf "%s\n" "$ROWS" | awk -F'|' '{print $4}' | while IFS= read -r p; do basename "$(dirname "$p")"; done | paste -sd',' -)"
      CMD+=(--sample-ids "$IDS")
    fi
  fi
  CMD+=(--block-length "$(prompt "Bootstrap block length (frames; test sensitivity)" "50")"
        --bootstrap "$(prompt "Bootstrap replicates (0=disabled, otherwise >=50)" "200")"
        --seed "$(prompt "Random seed" "0")")
fi
CMD+=(--progress)
printf 'Running: '; printf '%q ' "${CMD[@]}"; echo
"${CMD[@]}"
