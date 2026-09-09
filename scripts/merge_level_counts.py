#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median

CASE_DISTANCE_RE = re.compile(r"_d(?P<distance>\d+)_")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Summarize fixed level_params.csv into cluster-count distribution summaries.")
    p.add_argument("--fixed-level-params-csv", required=True, type=Path)
    p.add_argument("--summary-csv", required=True, type=Path)
    p.add_argument("--histogram-csv", required=True, type=Path)
    return p.parse_args()


def parse_distance(case_name: str) -> int:
    m = CASE_DISTANCE_RE.search(case_name)
    if not m:
        raise ValueError(f"could not parse distance from case name: {case_name}")
    return int(m.group("distance"))


def percentile(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    idx = (len(ordered) - 1) * q
    lo = math.floor(idx)
    hi = math.ceil(idx)
    if lo == hi:
        return ordered[lo]
    frac = idx - lo
    return ordered[lo] * (1.0 - frac) + ordered[hi] * frac


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    shot_keys_by_case_method: dict[tuple[str, str], set[tuple[int, int]]] = defaultdict(set)
    levels_by_case_method: dict[tuple[str, str], set[int]] = defaultdict(set)
    counts: dict[tuple[str, str, int, int, int], int] = {}
    assigned: dict[tuple[str, str, int, int, int], int] = {}
    max_active: dict[tuple[str, str, int, int, int], int] = {}
    max_events: dict[tuple[str, str, int, int, int], int] = {}
    max_critical: dict[tuple[str, str, int, int, int], int] = {}

    for row in rows:
        case = row["case"]
        method = row.get("method", "")
        level = int(row["level"])
        rep = int(row["rep"])
        shot = int(row["shot"])
        shot_key = (rep, shot)
        key = (case, method, level, rep, shot)
        cm = (case, method)
        shot_keys_by_case_method[cm].add(shot_key)
        levels_by_case_method[cm].add(level)
        counts[key] = int(float(row.get("cluster_count", "0") or 0))
        assigned[key] = int(float(row.get("assigned_active_detector_count", "0") or 0))
        max_active[key] = int(float(row.get("max_cluster_active_detectors", "0") or 0))
        max_events[key] = int(float(row.get("max_cluster_event_count", "0") or 0))
        max_critical[key] = int(float(row.get("max_cluster_critical_path_event_count", "0") or 0))

    summary_rows: list[dict[str, object]] = []
    hist_rows: list[dict[str, object]] = []
    for case, method in sorted(shot_keys_by_case_method.keys(), key=lambda x: (parse_distance(x[0]), x[1])):
        shot_keys = sorted(shot_keys_by_case_method[(case, method)])
        max_level = max(levels_by_case_method[(case, method)] or {0})
        for level in range(max_level + 1):
            cluster_values: list[float] = []
            assigned_values: list[float] = []
            max_active_values: list[float] = []
            max_event_values: list[float] = []
            max_critical_values: list[float] = []
            hist: Counter[int] = Counter()
            for rep, shot in shot_keys:
                key = (case, method, level, rep, shot)
                c = counts.get(key, 0)
                cluster_values.append(float(c))
                assigned_values.append(float(assigned.get(key, 0)))
                max_active_values.append(float(max_active.get(key, 0)))
                max_event_values.append(float(max_events.get(key, 0)))
                max_critical_values.append(float(max_critical.get(key, 0)))
                hist[c] += 1
            nshots = len(shot_keys)
            if nshots == 0:
                continue
            summary_rows.append({
                "case": case,
                "distance": parse_distance(case),
                "method": method,
                "level": level,
                "shots": nshots,
                "nonempty_shots": int(sum(1 for v in cluster_values if v > 0)),
                "nonempty_fraction": f"{sum(1 for v in cluster_values if v > 0) / nshots:.9f}",
                "mean_cluster_count": f"{mean(cluster_values):.9f}",
                "median_cluster_count": f"{median(cluster_values):.9f}",
                "p90_cluster_count": f"{percentile(cluster_values, 0.90):.9f}",
                "p95_cluster_count": f"{percentile(cluster_values, 0.95):.9f}",
                "p99_cluster_count": f"{percentile(cluster_values, 0.99):.9f}",
                "max_cluster_count": int(max(cluster_values)),
                "mean_assigned_active_detector_count": f"{mean(assigned_values):.9f}",
                "p95_assigned_active_detector_count": f"{percentile(assigned_values, 0.95):.9f}",
                "max_assigned_active_detector_count": int(max(assigned_values)),
                "mean_max_cluster_active_detectors": f"{mean(max_active_values):.9f}",
                "p95_max_cluster_active_detectors": f"{percentile(max_active_values, 0.95):.9f}",
                "max_cluster_active_detectors": int(max(max_active_values)),
                "mean_max_cluster_event_count": f"{mean(max_event_values):.9f}",
                "max_cluster_event_count": int(max(max_event_values)),
                "mean_max_cluster_critical_path_event_count": f"{mean(max_critical_values):.9f}",
                "max_cluster_critical_path_event_count": int(max(max_critical_values)),
            })
            for cluster_count in sorted(hist):
                shots = hist[cluster_count]
                hist_rows.append({
                    "case": case,
                    "distance": parse_distance(case),
                    "method": method,
                    "level": level,
                    "cluster_count": cluster_count,
                    "shots": shots,
                    "fraction": f"{shots / nshots:.9f}",
                })
    return summary_rows, hist_rows


def main() -> None:
    args = parse_args()
    rows = read_rows(args.fixed_level_params_csv)
    summary_rows, hist_rows = summarize(rows)
    write_csv(
        args.summary_csv,
        [
            "case", "distance", "method", "level", "shots", "nonempty_shots", "nonempty_fraction",
            "mean_cluster_count", "median_cluster_count", "p90_cluster_count", "p95_cluster_count",
            "p99_cluster_count", "max_cluster_count",
            "mean_assigned_active_detector_count", "p95_assigned_active_detector_count",
            "max_assigned_active_detector_count", "mean_max_cluster_active_detectors",
            "p95_max_cluster_active_detectors", "max_cluster_active_detectors",
            "mean_max_cluster_event_count", "max_cluster_event_count",
            "mean_max_cluster_critical_path_event_count", "max_cluster_critical_path_event_count",
        ],
        summary_rows,
    )
    write_csv(
        args.histogram_csv,
        ["case", "distance", "method", "level", "cluster_count", "shots", "fraction"],
        hist_rows,
    )


if __name__ == "__main__":
    main()
