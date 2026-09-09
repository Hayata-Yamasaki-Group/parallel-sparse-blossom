#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$ROOT"
BUILD="$ROOT/build/obj"
OUT="${OUT:-$ROOT/bin/parallel_eventcount_bench}"
CXX="${CXX:-g++}"
OPT="${OPT:--O3}"
mkdir -p "$BUILD" "$(dirname "$OUT")"

STAMP_FILE="$BUILD/.toolchain_stamp"
TOOLCHAIN_ID="$("$CXX" --version | head -n 1 || echo "$CXX unknown")"
PLATFORM_ID="$(uname -srm || echo unknown-platform)"
BUILD_SIGNATURE="$CXX|$TOOLCHAIN_ID|$PLATFORM_ID|$OPT|c++20"
if [[ ! -f "$STAMP_FILE" || "$(cat "$STAMP_FILE")" != "$BUILD_SIGNATURE" ]]; then
  echo "toolchain or platform changed; clearing stale objects in $BUILD"
  find "$BUILD" -maxdepth 1 -type f \( -name '*.o' -o -name '.toolchain_stamp' \) -delete
  printf '%s\n' "$BUILD_SIGNATURE" > "$STAMP_FILE"
fi

COMMON_FLAGS=(-std=c++20 "$OPT" -DNDEBUG -pthread -I"$SRC/src")
SOURCES=(
  src/pymatching/rand/rand_gen.cc
  src/pymatching/sparse_blossom/driver/mwpm_decoding.cc
  src/pymatching/sparse_blossom/driver/user_graph.cc
  src/pymatching/sparse_blossom/flooder/detector_node.cc
  src/pymatching/sparse_blossom/flooder/graph.cc
  src/pymatching/sparse_blossom/flooder/graph_fill_region.cc
  src/pymatching/sparse_blossom/flooder/graph_flooder.cc
  src/pymatching/sparse_blossom/flooder/match.cc
  src/pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.cc
  src/pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.cc
  src/pymatching/sparse_blossom/flooder_matcher_interop/region_edge.cc
  src/pymatching/sparse_blossom/matcher/alternating_tree.cc
  src/pymatching/sparse_blossom/matcher/mwpm.cc
  src/pymatching/sparse_blossom/parallel/growing_only_clusterer.cc
  src/pymatching/sparse_blossom/parallel/mwpm_live_snapshot.cc
  src/pymatching/sparse_blossom/parallel/mwpm_stopped_state.cc
  src/pymatching/sparse_blossom/parallel/processing_cluster.cc
  src/pymatching/sparse_blossom/parallel/cluster_subgraph.cc
  src/pymatching/sparse_blossom/parallel/lockstep_scheduler.cc
  src/pymatching/sparse_blossom/search/search_detector_node.cc
  src/pymatching/sparse_blossom/search/search_flooder.cc
  src/pymatching/sparse_blossom/search/search_graph.cc
)
OBJS=()
for i in "${!SOURCES[@]}"; do
  src_file="$SRC/${SOURCES[$i]}"
  obj="$BUILD/obj_$i.o"
  if [[ ! -e "$obj" || "$src_file" -nt "$obj" ]]; then
    echo "compile ${SOURCES[$i]}"
    "$CXX" "${COMMON_FLAGS[@]}" -c "$src_file" -o "$obj"
  fi
  OBJS+=("$obj")
done
bench_obj="$BUILD/bench.o"
if [[ ! -e "$bench_obj" || "$SRC/src/local_parallel_bench.cc" -nt "$bench_obj" ]]; then
  echo "compile src/local_parallel_bench.cc"
  "$CXX" "${COMMON_FLAGS[@]}" -c "$SRC/src/local_parallel_bench.cc" -o "$bench_obj"
fi
"$CXX" -pthread "${OBJS[@]}" "$bench_obj" -o "$OUT"
echo "Built $OUT"
