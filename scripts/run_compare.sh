#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/scripts/lib.sh"

P_VALUE="${P_VALUE:-0.0001}"
SHOTS="${SHOTS:-256}"
MAX_SHOTS="${MAX_SHOTS:-$SHOTS}"
DISTANCES="${DISTANCES:-29}"
WORKERS="${WORKERS:-1}"
REPS="${REPS:-1}"
SEED_BASE="${SEED_BASE:-200}"
NOISE_MODEL="${NOISE_MODEL:-physical}"
L_RATIO="${L_RATIO:-1.0}"
MIN_BUFFER_RATIO="${MIN_BUFFER_RATIO:-2.2}"
PHI_FLOOR="${PHI_FLOOR:-0.01}"
PHI_Q="${PHI_Q:-0.10}"
PARAMETER_SCHEDULE="${PARAMETER_SCHEDULE:-1}"
MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY=1
ALL_PAIRS_MAX_NODES="${ALL_PAIRS_MAX_NODES:-}"
ALL_PAIRS_CACHE_MAX_GIB="${ALL_PAIRS_CACHE_MAX_GIB:-}"
CASE_DIR="${CASE_DIR:-}"
CASE_TEMPLATE="${CASE_TEMPLATE:-}"
OUT_DIR="${OUT_DIR:-}"
GENERATE_CASES=1
AUTO_BUILD=1
FIXED_BOUNDS="${PYMATCHING_FIXED_BOUNDS:-${FIXED_BOUNDS:-}}"
ASSIGNMENT_LOG="${PYMATCHING_ASSIGNMENT_LOG:-${ASSIGNMENT_LOG:-0}}"
CLUSTER_LOG="${PYMATCHING_CLUSTER_LOG:-${CLUSTER_LOG:-0}}"
MISTAKE_LOG="${PYMATCHING_MISTAKE_LOG:-${MISTAKE_LOG:-0}}"
SWEEP_PLOT="${PYMATCHING_SWEEP_PLOT:-${SWEEP_PLOT:-0}}"

_FIXED_LABEL_WAS_SET=${FIXED_LABEL+x}
FIXED_LABEL="${FIXED_LABEL:-fixed ideal}"
if [[ -n "$FIXED_BOUNDS" && -z "${_FIXED_LABEL_WAS_SET:-}" ]]; then
  FIXED_LABEL="fixed calibrated"
fi
PLOT_TITLE="${PLOT_TITLE:-Global vs Fixed Event Counts}"

usage() {
  cat <<'EOF_USAGE'
usage: bash ./scripts/run_compare.sh [options]

Runs fixed clustering on the selected cases, then writes a CSV and SVG comparing:
  - global sparse blossom event count
  - fixed-clustering ideal event count

Common options:
  --case-dir DIR
  --case-template TEMPLATE
  --distances SPEC
  --workers N
  --shots N
  --max-shots N
  --reps N
  --seed-base N
  --p VALUE
  --noise-model physical|uniform
  --out-dir DIR
  --l-ratio VALUE
  --min-buffer-ratio VALUE
  --min-buffer-ratio-first-level-only
  --fixed-bounds D1:B1[,D2:B2,...]
  --assignment-log 0|1        write heavy cluster member / detector assignment logs
  --cluster-log 0|1           write per-cluster event logs
  --mistake-log 0|1           write mistake-only shot/cluster/detector logs
  --sweep-plot 0|1            keep per-method plot/debug outputs
  --phi-sequence-floor VALUE
  --phi-sequence-q VALUE
  --parameter-schedule       use paper Algorithm 3 schedule (default)
  --phi-sequence             use the legacy phi-sequence schedule instead
  --all-pairs-max-nodes N
  --all-pairs-cache-max-gib VALUE
  --no-generate-cases
  --no-build

Outputs:
  results/compare_<timestamp>/compare.csv
  results/compare_<timestamp>/compare.svg
  results/compare_<timestamp>/config.env
  results/compare_<timestamp>/level_cluster_count_summary.csv
  results/compare_<timestamp>/level_cluster_count_histogram.csv
EOF_USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case-dir) require_value "$1" "${2-}"; CASE_DIR="$2"; shift 2 ;;
    --case-template) require_value "$1" "${2-}"; CASE_TEMPLATE="$2"; shift 2 ;;
    --distances) require_value "$1" "${2-}"; DISTANCES="$2"; shift 2 ;;
    --workers) require_value "$1" "${2-}"; WORKERS="$2"; shift 2 ;;
    --shots) require_value "$1" "${2-}"; SHOTS="$2"; shift 2 ;;
    --max-shots) require_value "$1" "${2-}"; MAX_SHOTS="$2"; shift 2 ;;
    --reps) require_value "$1" "${2-}"; REPS="$2"; shift 2 ;;
    --seed-base) require_value "$1" "${2-}"; SEED_BASE="$2"; shift 2 ;;
    --p) require_value "$1" "${2-}"; P_VALUE="$2"; shift 2 ;;
    --noise-model)
      require_value "$1" "${2-}"
      NOISE_MODEL="$2"
      case "$NOISE_MODEL" in
        physical|physically_motivated) NOISE_MODEL="physical" ;;
        uniform) ;;
        *) die "--noise-model must be 'physical' or 'uniform'" ;;
      esac
      shift 2
      ;;
    --out-dir) require_value "$1" "${2-}"; OUT_DIR="$2"; shift 2 ;;
    --l-ratio) require_value "$1" "${2-}"; L_RATIO="$2"; shift 2 ;;
    --min-buffer-ratio) require_value "$1" "${2-}"; MIN_BUFFER_RATIO="$2"; shift 2 ;;
    --min-buffer-ratio-first-level-only) MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY=1; shift ;;
    --fixed-bounds) require_value "$1" "${2-}"; FIXED_BOUNDS="$2"; shift 2 ;;
    --assignment-log) require_value "$1" "${2-}"; ASSIGNMENT_LOG="$2"; shift 2 ;;
    --cluster-log) require_value "$1" "${2-}"; CLUSTER_LOG="$2"; shift 2 ;;
    --mistake-log) require_value "$1" "${2-}"; MISTAKE_LOG="$2"; shift 2 ;;
    --sweep-plot) require_value "$1" "${2-}"; SWEEP_PLOT="$2"; shift 2 ;;
    --phi-sequence-floor) require_value "$1" "${2-}"; PHI_FLOOR="$2"; shift 2 ;;
    --phi-sequence-q) require_value "$1" "${2-}"; PHI_Q="$2"; shift 2 ;;
    --parameter-schedule) PARAMETER_SCHEDULE=1; shift ;;
    --phi-sequence) PARAMETER_SCHEDULE=0; shift ;;
    --all-pairs-max-nodes) require_value "$1" "${2-}"; ALL_PAIRS_MAX_NODES="$2"; shift 2 ;;
    --all-pairs-cache-max-gib) require_value "$1" "${2-}"; ALL_PAIRS_CACHE_MAX_GIB="$2"; shift 2 ;;
    --no-generate-cases) GENERATE_CASES=0; shift ;;
    --no-build) AUTO_BUILD=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option '$1'" ;;
  esac
done

if [[ -z "$CASE_DIR" ]]; then
  CASE_DIR="cases/default_p1e4_d5_d${DISTANCES}_s${SHOTS}"
fi
if [[ -z "$CASE_TEMPLATE" ]]; then
  CASE_TEMPLATE="case_r{distance}_d{distance}_p{p4}_s${SHOTS}_seed{seed}_full_idle_with_brd.txt"
fi
[[ "$CASE_DIR" != /* ]] && CASE_DIR="$ROOT/$CASE_DIR"
if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="$ROOT/results/compare_$(date +%Y%m%d_%H%M%S)"
elif [[ "$OUT_DIR" != /* ]]; then
  OUT_DIR="$ROOT/$OUT_DIR"
fi
mkdir -p "$OUT_DIR"
FIXED_OUT="$OUT_DIR/fixed"
COMPARE_CSV="$OUT_DIR/compare.csv"
COMPARE_SVG="$OUT_DIR/compare.svg"
RUN_CONFIG="$OUT_DIR/config.env"
LEVEL_COUNT_SUMMARY_CSV="$OUT_DIR/level_cluster_count_summary.csv"
LEVEL_COUNT_HISTOGRAM_CSV="$OUT_DIR/level_cluster_count_histogram.csv"

if [[ "$GENERATE_CASES" -eq 1 ]]; then
  python3 -m pip show stim >/dev/null 2>&1 || python3 -m pip install stim
  bash "$ROOT/scripts/gen_cases.sh" --case-dir "$CASE_DIR" --distances "$DISTANCES" --p "$P_VALUE" --shots "$SHOTS" --seed-base "$SEED_BASE" --noise-model "$NOISE_MODEL"
fi
if [[ "$AUTO_BUILD" -eq 1 ]]; then
  bash "$ROOT/scripts/build.sh"
fi

common_args=(
  --case-dir "$CASE_DIR"
  --case-template "$CASE_TEMPLATE"
  --distances "$DISTANCES"
  --p "$P_VALUE"
  --shots "$SHOTS"
  --max-shots "$MAX_SHOTS"
  --workers "$WORKERS"
  --reps "$REPS"
  --seed-base "$SEED_BASE"
  --l-ratio "$L_RATIO"
  --min-buffer-ratio "$MIN_BUFFER_RATIO"
  --phi-sequence-floor "$PHI_FLOOR"
  --phi-sequence-q "$PHI_Q"
)
if [[ "$PARAMETER_SCHEDULE" -eq 1 ]]; then
  common_args+=(--parameter-schedule)
else
  common_args+=(--phi-sequence)
fi
[[ "$MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY" -eq 1 ]] && common_args+=(--min-buffer-ratio-first-level-only)
[[ -n "$ALL_PAIRS_MAX_NODES" ]] && common_args+=(--all-pairs-max-nodes "$ALL_PAIRS_MAX_NODES")
[[ -n "$ALL_PAIRS_CACHE_MAX_GIB" ]] && common_args+=(--all-pairs-cache-max-gib "$ALL_PAIRS_CACHE_MAX_GIB")
common_args+=(--assignment-log "$ASSIGNMENT_LOG")
common_args+=(--cluster-log "$CLUSTER_LOG")
common_args+=(--mistake-log "$MISTAKE_LOG")
common_args+=(--sweep-plot "$SWEEP_PLOT")
[[ "$AUTO_BUILD" -eq 0 ]] && common_args+=(--no-build)

cat >"$RUN_CONFIG" <<EOF_CONFIG
CASE_DIR=$CASE_DIR
CASE_TEMPLATE=$CASE_TEMPLATE
DISTANCES=$DISTANCES
P_VALUE=$P_VALUE
SHOTS=$SHOTS
MAX_SHOTS=$MAX_SHOTS
WORKERS=$WORKERS
REPS=$REPS
SEED_BASE=$SEED_BASE
NOISE_MODEL=$NOISE_MODEL
L_RATIO=$L_RATIO
MIN_BUFFER_RATIO=$MIN_BUFFER_RATIO
PHI_FLOOR=$PHI_FLOOR
PHI_Q=$PHI_Q
PARAMETER_SCHEDULE=$PARAMETER_SCHEDULE
ALL_PAIRS_MAX_NODES=$ALL_PAIRS_MAX_NODES
ALL_PAIRS_CACHE_MAX_GIB=$ALL_PAIRS_CACHE_MAX_GIB
FIXED_BOUNDS=${FIXED_BOUNDS:-}
ASSIGNMENT_LOG=$ASSIGNMENT_LOG
CLUSTER_LOG=$CLUSTER_LOG
MISTAKE_LOG=$MISTAKE_LOG
SWEEP_PLOT=$SWEEP_PLOT
EOF_CONFIG

# Fixed clustering.
fixed_args=(
  "${common_args[@]}"
  --cluster-from-match-paths 1
  --out-dir "$FIXED_OUT"
)
[[ -n "$FIXED_BOUNDS" ]] && fixed_args+=(--cluster-bounds "$FIXED_BOUNDS")
env \
  PYMATCHING_EVENTCOUNT_METHOD_LABEL=fixed \
  PYMATCHING_ADAPTIVE_CLUSTER_BOUNDS=0 \
  PYMATCHING_ADAPTIVE_BREAKPOINT_BOUNDS=0 \
  bash "$ROOT/scripts/run_sweep.sh" \
    "${fixed_args[@]}"

python3 "$ROOT/scripts/plot_compare.py" \
  --fixed-summary-csv "$FIXED_OUT/eventcount_summary.csv" \
  --comparison-csv "$COMPARE_CSV" \
  --plot-svg "$COMPARE_SVG" \
  --fixed-label "$FIXED_LABEL" \
  --plot-title "$PLOT_TITLE"
python3 "$ROOT/scripts/merge_level_counts.py" \
  --fixed-level-params-csv "$FIXED_OUT/level_params.csv" \
  --summary-csv "$LEVEL_COUNT_SUMMARY_CSV" \
  --histogram-csv "$LEVEL_COUNT_HISTOGRAM_CSV"

echo "Wrote $OUT_DIR"
echo "Wrote $COMPARE_CSV"
echo "Wrote $COMPARE_SVG"
echo "Wrote $RUN_CONFIG"
echo "Wrote $LEVEL_COUNT_SUMMARY_CSV"
echo "Wrote $LEVEL_COUNT_HISTOGRAM_CSV"
