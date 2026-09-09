#!/usr/bin/env bash

die() {
  echo "error: $*" >&2
  exit 1
}

require_value() {
  local flag="$1"
  local value="${2-}"
  [[ -n "$value" ]] || die "$flag requires a value"
}

parse_distances() {
  local spec="$1"
  local default_start="${2:-5}"
  local -a values=()
  if [[ "$spec" == *","* ]]; then
    IFS=',' read -r -a values <<<"$spec"
    for value in "${values[@]}"; do
      [[ "$value" =~ ^[0-9]+$ ]] || die "invalid distance '$value' in list '$spec'"
      echo "$value"
    done
    return 0
  fi
  if [[ "$spec" =~ ^([0-9]+):([0-9]+)(:([0-9]+))?$ ]]; then
    local start="${BASH_REMATCH[1]}"
    local stop="${BASH_REMATCH[2]}"
    local step="${BASH_REMATCH[4]:-2}"
    (( step > 0 )) || die "distance step must be positive"
    for ((d = start; d <= stop; d += step)); do
      echo "$d"
    done
    return 0
  fi
  if [[ "$spec" =~ ^[0-9]+$ ]]; then
    local stop="$spec"
    if (( stop <= default_start )); then
      echo "$stop"
      return 0
    fi
    for ((d = default_start; d <= stop; d += 2)); do
      echo "$d"
    done
    return 0
  fi
  die "unsupported distance specification '$spec'"
}

format_probability() {
  local p="$1"
  printf "%.4f" "$p"
}

render_case_template() {
  local template="$1"
  local distance="$2"
  local rounds="$3"
  local seed="$4"
  local shots="$5"
  local p_raw="${6-}"
  local rendered="$template"
  local p_fixed=""
  if [[ -n "$p_raw" ]]; then
    p_fixed="$(format_probability "$p_raw")"
  fi
  rendered="${rendered//\{distance\}/$distance}"
  rendered="${rendered//\{rounds\}/$rounds}"
  rendered="${rendered//\{seed\}/$seed}"
  rendered="${rendered//\{shots\}/$shots}"
  if [[ -n "$p_raw" ]]; then
    rendered="${rendered//\{p\}/$p_raw}"
    rendered="${rendered//\{p4\}/$p_fixed}"
  fi
  echo "$rendered"
}
