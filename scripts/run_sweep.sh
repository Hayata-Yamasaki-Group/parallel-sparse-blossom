#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/scripts/lib.sh"

BIN="${BIN:-$ROOT/bin/parallel_eventcount_bench}"
CASE_DIR=""
CASE_TEMPLATE=""
DISTANCE_SPEC=""
OUT_DIR="${OUT_DIR:-$ROOT/results/sweep_$(date +%Y%m%d_%H%M%S)}"
WORKERS=1
REPS=1
MAX_SHOTS=""
SEED_BASE=200
SHOTS=256
P_VALUE=""
DEFAULT_START_DISTANCE=5
L_RATIO=1
MIN_BUFFER_RATIO=2.2
MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY=0
PHI_SEQUENCE=0
PHI_SEQUENCE_FLOOR="0.01"
PHI_SEQUENCE_Q="0.10"
PARAMETER_SCHEDULE=1
ALL_PAIRS_MAX_NODES=""
ALL_PAIRS_CACHE_MAX_GIB=""
CLUSTER_FROM_MATCH_PATHS=1
CLUSTER_BOUNDS="${PYMATCHING_CLUSTER_BOUNDS:-}"
ASSIGNMENT_LOG="${PYMATCHING_ASSIGNMENT_LOG:-0}"
CLUSTER_LOG="${PYMATCHING_CLUSTER_LOG:-0}"
MISTAKE_LOG="${PYMATCHING_MISTAKE_LOG:-1}"
LEVEL_PARAMS_LOG="${PYMATCHING_LEVEL_PARAMS_LOG:-1}"
SWEEP_PLOT="${PYMATCHING_SWEEP_PLOT:-0}"
AUTO_BUILD=1

usage() {
  cat <<'EOF_USAGE'
usage: bash ./scripts/run_sweep.sh --case-dir DIR --case-template TEMPLATE --distances SPEC [options]

This is the internal event-count sweep runner used by run_compare.sh.

Common options:
  --out-dir DIR
  --workers N
  --reps N
  --max-shots N
  --seed-base N
  --shots N
  --p VALUE
  --l-ratio VALUE
  --min-buffer-ratio VALUE
  --min-buffer-ratio-first-level-only
  --phi-sequence
  --phi-sequence-floor VALUE
  --phi-sequence-q VALUE
  --parameter-schedule       use paper Algorithm 3 schedule (default)
  --all-pairs-max-nodes N
  --all-pairs-cache-max-gib VALUE
  --cluster-from-match-paths 0|1
  --cluster-bounds D1:B1[,D2:B2,...]
  --assignment-log 0|1        write cluster_members.csv and detector_assignments.csv
  --cluster-log 0|1           write cluster_events.csv and level_cluster_event_stats.csv
  --level-params-log 0|1      write level_params.csv used for core-count distributions
  --sweep-plot 0|1           keep per-method eventcount_plot_data/svg/fit outputs
  --no-build
EOF_USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case-dir) require_value "$1" "${2-}"; CASE_DIR="$2"; shift 2 ;;
    --case-template) require_value "$1" "${2-}"; CASE_TEMPLATE="$2"; shift 2 ;;
    --distances) require_value "$1" "${2-}"; DISTANCE_SPEC="$2"; shift 2 ;;
    --out-dir) require_value "$1" "${2-}"; OUT_DIR="$2"; shift 2 ;;
    --workers) require_value "$1" "${2-}"; WORKERS="$2"; shift 2 ;;
    --reps) require_value "$1" "${2-}"; REPS="$2"; shift 2 ;;
    --max-shots) require_value "$1" "${2-}"; MAX_SHOTS="$2"; shift 2 ;;
    --seed-base) require_value "$1" "${2-}"; SEED_BASE="$2"; shift 2 ;;
    --shots) require_value "$1" "${2-}"; SHOTS="$2"; shift 2 ;;
    --p) require_value "$1" "${2-}"; P_VALUE="$2"; shift 2 ;;
    --distance-start) require_value "$1" "${2-}"; DEFAULT_START_DISTANCE="$2"; shift 2 ;;
    --l-ratio) require_value "$1" "${2-}"; L_RATIO="$2"; shift 2 ;;
    --min-buffer-ratio) require_value "$1" "${2-}"; MIN_BUFFER_RATIO="$2"; shift 2 ;;
    --min-buffer-ratio-first-level-only) MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY=1; shift ;;
    --phi-sequence) PHI_SEQUENCE=1; PARAMETER_SCHEDULE=0; shift ;;
    --phi-sequence-floor) require_value "$1" "${2-}"; PHI_SEQUENCE_FLOOR="$2"; shift 2 ;;
    --phi-sequence-q) require_value "$1" "${2-}"; PHI_SEQUENCE_Q="$2"; shift 2 ;;
    --parameter-schedule) PARAMETER_SCHEDULE=1; PHI_SEQUENCE=0; shift ;;
    --all-pairs-max-nodes) require_value "$1" "${2-}"; ALL_PAIRS_MAX_NODES="$2"; shift 2 ;;
    --all-pairs-cache-max-gib) require_value "$1" "${2-}"; ALL_PAIRS_CACHE_MAX_GIB="$2"; shift 2 ;;
    --cluster-from-match-paths) require_value "$1" "${2-}"; CLUSTER_FROM_MATCH_PATHS="$2"; shift 2 ;;
    --cluster-bounds) require_value "$1" "${2-}"; CLUSTER_BOUNDS="$2"; shift 2 ;;
    --assignment-log) require_value "$1" "${2-}"; ASSIGNMENT_LOG="$2"; shift 2 ;;
    --cluster-log) require_value "$1" "${2-}"; CLUSTER_LOG="$2"; shift 2 ;;
    --mistake-log) require_value "$1" "${2-}"; MISTAKE_LOG="$2"; shift 2 ;;
    --level-params-log) require_value "$1" "${2-}"; LEVEL_PARAMS_LOG="$2"; shift 2 ;;
    --sweep-plot) require_value "$1" "${2-}"; SWEEP_PLOT="$2"; shift 2 ;;
    --no-build) AUTO_BUILD=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option '$1'" ;;
  esac
done

[[ -n "$CASE_DIR" ]] || die "--case-dir is required"
[[ -n "$CASE_TEMPLATE" ]] || die "--case-template is required"
[[ -n "$DISTANCE_SPEC" ]] || die "--distances is required"
if [[ "$CASE_TEMPLATE" == *"{p}"* || "$CASE_TEMPLATE" == *"{p4}"* ]]; then
  [[ -n "$P_VALUE" ]] || die "--p is required when case-template uses {p} or {p4}"
fi

if [[ ! -x "$BIN" ]]; then
  [[ "$AUTO_BUILD" -eq 1 ]] || die "binary not found or not executable: $BIN"
  OUT="$BIN" "$ROOT/scripts/build.sh"
fi

mkdir -p "$OUT_DIR"
RAW_CSV="$OUT_DIR/bench_rows.csv"
RAW_LOG="$OUT_DIR/bench_stderr.log"
CLUSTER_CSV="$OUT_DIR/cluster_events.csv"
CLUSTER_MEMBERS_CSV="$OUT_DIR/cluster_members.csv"
DETECTOR_ASSIGNMENTS_CSV="$OUT_DIR/detector_assignments.csv"
MISTAKE_SHOTS_CSV="$OUT_DIR/mistake_shots.csv"
MISTAKE_CLUSTERS_CSV="$OUT_DIR/mistake_clusters.csv"
MISTAKE_DETECTORS_CSV="$OUT_DIR/mistake_detectors.csv"
MISTAKE_MATCHING_CSV="$OUT_DIR/mistake_matching.csv"
LEVEL_PARAMS_CSV="$OUT_DIR/level_params.csv"
SUMMARY_CSV="$OUT_DIR/eventcount_summary.csv"
LEVEL_STATS_CSV="$OUT_DIR/level_cluster_event_stats.csv"
LEVEL_COUNT_SUMMARY_CSV="$OUT_DIR/level_cluster_count_summary.csv"
LEVEL_COUNT_HISTOGRAM_CSV="$OUT_DIR/level_cluster_count_histogram.csv"
PLOT_DATA_CSV="$OUT_DIR/eventcount_plot_data.csv"
PLOT_SVG="$OUT_DIR/eventcount_distance_loglog.svg"
FIT_CSV="$OUT_DIR/eventcount_scaling_fit.csv"
CONFIG_ENV="$OUT_DIR/run_config.env"
rm -f "$RAW_CSV" "$RAW_LOG" "$CLUSTER_CSV" "$CLUSTER_MEMBERS_CSV" "$DETECTOR_ASSIGNMENTS_CSV" "$MISTAKE_SHOTS_CSV" "$MISTAKE_CLUSTERS_CSV" "$MISTAKE_DETECTORS_CSV" "$MISTAKE_MATCHING_CSV" "$LEVEL_PARAMS_CSV" "$SUMMARY_CSV" "$LEVEL_STATS_CSV" "$LEVEL_COUNT_SUMMARY_CSV" "$LEVEL_COUNT_HISTOGRAM_CSV" "$PLOT_DATA_CSV" "$PLOT_SVG" "$FIT_CSV" "$CONFIG_ENV"

cat >"$CONFIG_ENV" <<EOF_CONFIG
BIN=$BIN
CASE_DIR=$CASE_DIR
CASE_TEMPLATE=$CASE_TEMPLATE
DISTANCES=$DISTANCE_SPEC
WORKERS=$WORKERS
REPS=$REPS
MAX_SHOTS=${MAX_SHOTS:-}
SEED_BASE=$SEED_BASE
SHOTS=$SHOTS
PYMATCHING_CLUSTER_FROM_MATCH_PATHS=$CLUSTER_FROM_MATCH_PATHS
PYMATCHING_CLUSTER_BOUNDS=${CLUSTER_BOUNDS:-}
PYMATCHING_ASSIGNMENT_LOG=$ASSIGNMENT_LOG
PYMATCHING_CLUSTER_LOG=$CLUSTER_LOG
PYMATCHING_MISTAKE_LOG=$MISTAKE_LOG
PYMATCHING_LEVEL_PARAMS_LOG=$LEVEL_PARAMS_LOG
PYMATCHING_SWEEP_PLOT=$SWEEP_PLOT
PYMATCHING_CLUSTER_PATH_L_RATIO=$L_RATIO
PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO=$MIN_BUFFER_RATIO
PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY=$MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY
PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE=$PHI_SEQUENCE
PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_FLOOR=${PHI_SEQUENCE_FLOOR:-}
PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_Q=${PHI_SEQUENCE_Q:-}
PYMATCHING_CLUSTER_PATH_PARAMETER_SCHEDULE=$PARAMETER_SCHEDULE
PYMATCHING_PROCESSING_CLUSTER_ALL_PAIRS_MAX_NODES=${ALL_PAIRS_MAX_NODES:-}
PYMATCHING_PROCESSING_CLUSTER_CACHE_MAX_GIB=${ALL_PAIRS_CACHE_MAX_GIB:-}
EOF_CONFIG

mapfile -t DISTANCES < <(parse_distances "$DISTANCE_SPEC" "$DEFAULT_START_DISTANCE")
[[ "${#DISTANCES[@]}" -gt 0 ]] || die "no distances selected"
first_case=1
skipped_no_matching_paths=0
for distance in "${DISTANCES[@]}"; do
  rounds="$distance"
  seed=$((SEED_BASE + distance))
  case_name="$(render_case_template "$CASE_TEMPLATE" "$distance" "$rounds" "$seed" "$SHOTS" "$P_VALUE")"
  case_path="$CASE_DIR/$case_name"
  [[ -f "$case_path" ]] || die "case file not found: $case_path"
  echo "Running distance $distance with case $case_path"
  tmp_csv="$(mktemp)"
  tmp_log="$(mktemp)"
  env_args=(
    PYMATCHING_PARALLEL_WORKERS="$WORKERS"
    PYMATCHING_EVENT_COUNT_ONLY=1
    PYMATCHING_CLUSTER_FROM_MATCH_PATHS="$CLUSTER_FROM_MATCH_PATHS"
    PYMATCHING_CLUSTER_PATH_L_RATIO="$L_RATIO"
    PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO="$MIN_BUFFER_RATIO"
    PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY="$MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY"
    PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE="$PHI_SEQUENCE"
    PYMATCHING_CLUSTER_PATH_PARAMETER_SCHEDULE="$PARAMETER_SCHEDULE"
    PYMATCHING_COMPONENT_IDS_ONLY_GRAPH_CACHE=1
    PYMATCHING_SKIP_RADIUS_NEIGHBOR_PRECOMPUTE=1
    PYMATCHING_SKIP_FORCED_MAX_EXACT_CLUSTER_METADATA=1
    PYMATCHING_DISABLE_WALL_TIMING=1
    PYMATCHING_EVENTCOUNT_ONLY_NO_WALL_TIMING=1
  )

  if [[ "$CLUSTER_LOG" == "1" ]]; then
    env_args+=(PYMATCHING_EVENTCOUNT_CLUSTER_CSV="$CLUSTER_CSV")
  fi
  if [[ "$ASSIGNMENT_LOG" == "1" ]]; then
    env_args+=(PYMATCHING_CLUSTER_MEMBERS_CSV="$CLUSTER_MEMBERS_CSV")
    env_args+=(PYMATCHING_DETECTOR_ASSIGNMENTS_CSV="$DETECTOR_ASSIGNMENTS_CSV")
  fi
  if [[ "$MISTAKE_LOG" == "1" ]]; then
    env_args+=(PYMATCHING_MISTAKE_SHOTS_CSV="$MISTAKE_SHOTS_CSV")
    env_args+=(PYMATCHING_MISTAKE_CLUSTERS_CSV="$MISTAKE_CLUSTERS_CSV")
    env_args+=(PYMATCHING_MISTAKE_DETECTORS_CSV="$MISTAKE_DETECTORS_CSV")
    env_args+=(PYMATCHING_MISTAKE_MATCHING_CSV="$MISTAKE_MATCHING_CSV")
  fi
  if [[ "$LEVEL_PARAMS_LOG" == "1" ]]; then
    env_args+=(PYMATCHING_LEVEL_PARAMS_CSV="$LEVEL_PARAMS_CSV")
  fi
  [[ -n "$MAX_SHOTS" ]] && env_args+=(PYMATCHING_BENCH_MAX_SHOTS="$MAX_SHOTS")
  [[ -n "$PHI_SEQUENCE_FLOOR" ]] && env_args+=(PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_FLOOR="$PHI_SEQUENCE_FLOOR")
  [[ -n "$PHI_SEQUENCE_Q" ]] && env_args+=(PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_Q="$PHI_SEQUENCE_Q")
  [[ -n "$ALL_PAIRS_MAX_NODES" ]] && env_args+=(PYMATCHING_PROCESSING_CLUSTER_ALL_PAIRS_MAX_NODES="$ALL_PAIRS_MAX_NODES")
  [[ -n "$ALL_PAIRS_CACHE_MAX_GIB" ]] && env_args+=(PYMATCHING_PROCESSING_CLUSTER_CACHE_MAX_GIB="$ALL_PAIRS_CACHE_MAX_GIB")
  [[ -n "$CLUSTER_BOUNDS" ]] && env_args+=(PYMATCHING_CLUSTER_BOUNDS="$CLUSTER_BOUNDS")

  run_env=(env)
  if [[ -z "$CLUSTER_BOUNDS" ]]; then
    # Avoid accidentally inheriting an outer debug bounds override.
    run_env+=(-u PYMATCHING_CLUSTER_BOUNDS)
  fi
  if "${run_env[@]}" "${env_args[@]}" "$BIN" "$REPS" "$case_path" >"$tmp_csv" 2>"$tmp_log"; then
    :
  else
    rc=$?
    if grep -qF \
      "Could not infer a path-distance unit because no matching paths were measured." \
      "$tmp_log"; then
      {
        echo "### distance=$distance case=$case_path SKIPPED reason=no_matching_paths"
        cat "$tmp_log"
        echo
      } >>"$RAW_LOG"
      echo "Skipping distance $distance: no matching paths were measured." >&2
      skipped_no_matching_paths=$((skipped_no_matching_paths + 1))
      rm -f "$tmp_csv" "$tmp_log"
      continue
    fi

    {
      echo "### distance=$distance case=$case_path FAILED rc=$rc"
      cat "$tmp_log"
      echo
    } >>"$RAW_LOG"
    cat "$tmp_log" >&2
    rm -f "$tmp_csv" "$tmp_log"
    exit "$rc"
  fi
  if [[ "$first_case" -eq 1 ]]; then
    cp "$tmp_csv" "$RAW_CSV"
    first_case=0
  else
    tail -n +2 "$tmp_csv" >>"$RAW_CSV"
  fi
  { echo "### distance=$distance case=$case_path"; cat "$tmp_log"; echo; } >>"$RAW_LOG"
  rm -f "$tmp_csv" "$tmp_log"
done

if [[ "$first_case" -eq 1 ]]; then
  die "all selected distances were skipped because no matching paths were measured"
fi
if [[ "$skipped_no_matching_paths" -gt 0 ]]; then
  echo "Skipped $skipped_no_matching_paths distance(s) with no measured matching paths." >&2
fi

python3 "$ROOT/scripts/analyze.py" \
  --bench-csv "$RAW_CSV" \
  --cluster-csv "$CLUSTER_CSV" \
  --level-params-csv "$LEVEL_PARAMS_CSV" \
  --summary-csv "$SUMMARY_CSV" \
  --level-stats-csv "$LEVEL_STATS_CSV" \
  --level-count-summary-csv "$LEVEL_COUNT_SUMMARY_CSV" \
  --level-count-histogram-csv "$LEVEL_COUNT_HISTOGRAM_CSV" \
  --plot-data-csv "$PLOT_DATA_CSV" \
  --plot-svg "$PLOT_SVG" \
  --fit-csv "$FIT_CSV"

if [[ "$CLUSTER_LOG" != "1" ]]; then
  rm -f "$CLUSTER_CSV" "$LEVEL_STATS_CSV"
fi
if [[ "$SWEEP_PLOT" != "1" ]]; then
  rm -f "$PLOT_DATA_CSV" "$PLOT_SVG" "$FIT_CSV"
fi

echo "Wrote $RAW_CSV"
echo "Wrote $RAW_LOG"
[[ "$CLUSTER_LOG" == "1" ]] && echo "Wrote $CLUSTER_CSV"
[[ "$ASSIGNMENT_LOG" == "1" ]] && echo "Wrote $CLUSTER_MEMBERS_CSV"
[[ "$ASSIGNMENT_LOG" == "1" ]] && echo "Wrote $DETECTOR_ASSIGNMENTS_CSV"
[[ "$MISTAKE_LOG" == "1" && -f "$MISTAKE_SHOTS_CSV" ]] && echo "Wrote $MISTAKE_SHOTS_CSV"
[[ "$MISTAKE_LOG" == "1" && -f "$MISTAKE_CLUSTERS_CSV" ]] && echo "Wrote $MISTAKE_CLUSTERS_CSV"
[[ "$MISTAKE_LOG" == "1" && -f "$MISTAKE_DETECTORS_CSV" ]] && echo "Wrote $MISTAKE_DETECTORS_CSV"
[[ "$MISTAKE_LOG" == "1" && -f "$MISTAKE_MATCHING_CSV" ]] && echo "Wrote $MISTAKE_MATCHING_CSV"
[[ "$LEVEL_PARAMS_LOG" == "1" ]] && echo "Wrote $LEVEL_PARAMS_CSV"
[[ "$LEVEL_PARAMS_LOG" == "1" ]] && echo "Wrote $LEVEL_COUNT_SUMMARY_CSV"
[[ "$LEVEL_PARAMS_LOG" == "1" ]] && echo "Wrote $LEVEL_COUNT_HISTOGRAM_CSV"
echo "Wrote $SUMMARY_CSV"
[[ -f "$PLOT_SVG" ]] && echo "Wrote $PLOT_SVG"

true
