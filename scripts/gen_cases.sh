#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "$ROOT/scripts/lib.sh"

GENERATOR="$ROOT/tools/generate_stim_case_with_idle.py"
CASE_DIR=""
DISTANCE_SPEC=""
P_VALUE=""
# Paper noise models.
NOISE_MODEL="physical"
P1Q_FACTOR="0.1"
IDLE_FACTOR="0"
SHOTS=256
SEED_BASE=200
DEFAULT_START_DISTANCE=5
ROUNDS_MODE="distance"
ROUNDS_VALUE=""
OUT_TEMPLATE='case_r{rounds}_d{distance}_p{p4}_s{shots}_seed{seed}_full_idle_with_brd.txt'

usage() {
  cat <<'EOF_USAGE'
usage: bash ./scripts/gen_cases.sh --case-dir DIR --distances SPEC --p VALUE [options]

Options:
  --shots N
  --seed-base N
  --noise-model physical|uniform
  --rounds N
  --distance-start N
  --out-template TEMPLATE
EOF_USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case-dir) require_value "$1" "${2-}"; CASE_DIR="$2"; shift 2 ;;
    --distances) require_value "$1" "${2-}"; DISTANCE_SPEC="$2"; shift 2 ;;
    --p) require_value "$1" "${2-}"; P_VALUE="$2"; shift 2 ;;
    --shots) require_value "$1" "${2-}"; SHOTS="$2"; shift 2 ;;
    --seed-base) require_value "$1" "${2-}"; SEED_BASE="$2"; shift 2 ;;
    --noise-model)
      require_value "$1" "${2-}"
      NOISE_MODEL="$2"
      case "$NOISE_MODEL" in
        physical|physically_motivated) P1Q_FACTOR="0.1"; IDLE_FACTOR="0" ;;
        uniform) P1Q_FACTOR="1.0"; IDLE_FACTOR="1.0" ;;
        *) die "--noise-model must be 'physical' or 'uniform'" ;;
      esac
      shift 2
      ;;
    --rounds) require_value "$1" "${2-}"; ROUNDS_MODE="fixed"; ROUNDS_VALUE="$2"; shift 2 ;;
    --distance-start) require_value "$1" "${2-}"; DEFAULT_START_DISTANCE="$2"; shift 2 ;;
    --out-template) require_value "$1" "${2-}"; OUT_TEMPLATE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option '$1'" ;;
  esac
done

[[ -n "$CASE_DIR" ]] || die "--case-dir is required"
[[ -n "$DISTANCE_SPEC" ]] || die "--distances is required"
[[ -n "$P_VALUE" ]] || die "--p is required"
[[ -f "$GENERATOR" ]] || die "generator not found: $GENERATOR"
mkdir -p "$CASE_DIR"
mapfile -t DISTANCES < <(parse_distances "$DISTANCE_SPEC" "$DEFAULT_START_DISTANCE")

for distance in "${DISTANCES[@]}"; do
  rounds="$distance"
  [[ "$ROUNDS_MODE" == "fixed" ]] && rounds="$ROUNDS_VALUE"
  seed=$((SEED_BASE + distance))
  out_name="$(render_case_template "$OUT_TEMPLATE" "$distance" "$rounds" "$seed" "$SHOTS" "$P_VALUE")"
  out_path="$CASE_DIR/$out_name"
  echo "Generating d=$distance -> $out_path"
  python3 "$GENERATOR" --distance "$distance" --rounds "$rounds" --p "$P_VALUE" --p1q-factor "$P1Q_FACTOR" --idle-factor "$IDLE_FACTOR" --shots "$SHOTS" --seed "$seed" --out "$out_path"
done
