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
    parser = argparse.ArgumentParser(description="Summarize event-count distance sweeps.")
    parser.add_argument("--bench-csv", required=True, type=Path)
    parser.add_argument("--cluster-csv", required=True, type=Path)
    parser.add_argument("--level-params-csv", type=Path)
    parser.add_argument("--summary-csv", required=True, type=Path)
    parser.add_argument("--level-stats-csv", required=True, type=Path)
    parser.add_argument("--level-count-summary-csv", type=Path)
    parser.add_argument("--level-count-histogram-csv", type=Path)
    parser.add_argument("--plot-data-csv", required=True, type=Path)
    parser.add_argument("--plot-svg", required=True, type=Path)
    parser.add_argument("--fit-csv", type=Path)
    parser.add_argument(
        "--plot-title",
        default="Global vs Parallel Sparse Blossom Events per Round",
    )
    return parser.parse_args()


def parse_distance(case_name: str) -> int:
    match = CASE_DISTANCE_RE.search(case_name)
    if not match:
        raise ValueError(f"could not parse distance from case name: {case_name}")
    return int(match.group("distance"))


def to_float(row: dict[str, str], key: str) -> float:
    raw = row.get(key, "")
    return float(raw) if raw not in ("", None) else 0.0


def first_float(row: dict[str, str], keys: list[str]) -> float:
    for key in keys:
        raw = row.get(key, "")
        if raw not in ("", None):
            return float(raw)
    return 0.0


def percentile(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = (len(ordered) - 1) * q
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize_bench_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    grouped: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        grouped[row["case"]][row["mode"]].append(row)

    summary_rows: list[dict[str, object]] = []
    for case_name, mode_rows in sorted(grouped.items(), key=lambda item: parse_distance(item[0])):
        global_rows = mode_rows.get("global", [])
        parallel_rows = [row for mode, rows_for_mode in mode_rows.items() if mode.startswith("parallel_") for row in rows_for_mode]
        if not global_rows or not parallel_rows:
            continue
        global_row = global_rows[0]
        parallel_row = parallel_rows[0]
        global_events = to_float(global_row, "global_mwpm_events_per_shot")
        global_events_sem = to_float(global_row, "global_mwpm_events_std_error_per_shot")
        parallel_ideal_events = first_float(
            parallel_row,
            ["parallel_level_critical_events_per_shot", "parallel_max_cluster_events_per_shot"],
        )
        parallel_ideal_events_sem = first_float(
            parallel_row,
            ["parallel_level_critical_events_std_error_per_shot", "parallel_max_cluster_events_std_error_per_shot"],
        )
        parallel_total_events = to_float(parallel_row, "parallel_cluster_events_total_per_shot")
        parallel_total_events_sem = to_float(parallel_row, "parallel_cluster_events_total_std_error_per_shot")
        parallel_algorithmic_time_sem = to_float(parallel_row, "parallel_critical_path_algorithmic_time_std_error_per_shot")
        speedup = global_events / parallel_ideal_events if parallel_ideal_events > 0 else 0.0
        distance = parse_distance(case_name)
        global_events_per_round = global_events / distance if distance > 0 else 0.0
        global_events_per_round_sem = global_events_sem / distance if distance > 0 else 0.0
        parallel_ideal_events_per_round = parallel_ideal_events / distance if distance > 0 else 0.0
        parallel_ideal_events_per_round_sem = parallel_ideal_events_sem / distance if distance > 0 else 0.0
        parallel_total_events_per_round = parallel_total_events / distance if distance > 0 else 0.0
        parallel_total_events_per_round_sem = parallel_total_events_sem / distance if distance > 0 else 0.0
        summary_rows.append(
            {
                "case": case_name,
                "distance": distance,
                "parallel_mode": parallel_row["mode"],
                "global_events_per_shot": f"{global_events:.6f}",
                "global_events_std_error_per_shot": f"{global_events_sem:.6f}",
                "parallel_ideal_events_per_shot": f"{parallel_ideal_events:.6f}",
                "parallel_ideal_events_std_error_per_shot": f"{parallel_ideal_events_sem:.6f}",
                "parallel_total_cluster_events_per_shot": f"{parallel_total_events:.6f}",
                "parallel_total_cluster_events_std_error_per_shot": f"{parallel_total_events_sem:.6f}",
                "global_events_per_shot_per_round": f"{global_events_per_round:.9f}",
                "global_events_per_shot_per_round_std_error": f"{global_events_per_round_sem:.9f}",
                "parallel_ideal_events_per_shot_per_round": f"{parallel_ideal_events_per_round:.9f}",
                "parallel_ideal_events_per_shot_per_round_std_error": f"{parallel_ideal_events_per_round_sem:.9f}",
                "parallel_total_cluster_events_per_shot_per_round": f"{parallel_total_events_per_round:.9f}",
                "parallel_total_cluster_events_per_shot_per_round_std_error": f"{parallel_total_events_per_round_sem:.9f}",
                "parallel_nonempty_event_levels_per_shot": f"{to_float(parallel_row, 'parallel_nonempty_event_levels_per_shot'):.6f}",
                "parallel_critical_path_algorithmic_time_per_shot": f"{to_float(parallel_row, 'parallel_critical_path_algorithmic_time_per_shot'):.6f}",
                "parallel_critical_path_algorithmic_time_std_error_per_shot": f"{parallel_algorithmic_time_sem:.6f}",
                "avg_clusters": f"{to_float(parallel_row, 'avg_clusters'):.6f}",
                "avg_max_level": f"{to_float(parallel_row, 'avg_max_level'):.6f}",
                "ideal_event_count_speedup_global_over_parallel": f"{speedup:.6f}",
                "mistakes": parallel_row.get("mistakes", ""),
                "exceptions": parallel_row.get("exceptions", ""),
            }
        )
    return summary_rows


def summarize_cluster_rows(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    shots_per_case: dict[str, set[tuple[int, int]]] = defaultdict(set)
    events_by_case_level: dict[tuple[str, int], list[float]] = defaultdict(list)
    active_by_case_level: dict[tuple[str, int], list[float]] = defaultdict(list)
    diameter_by_case_level: dict[tuple[str, int], list[float]] = defaultdict(list)
    buffer_by_case_level: dict[tuple[str, int], list[float]] = defaultdict(list)
    clusters_per_shot: dict[tuple[str, int], Counter[tuple[int, int]]] = defaultdict(Counter)

    for row in rows:
        case_name = row["case"]
        level = int(row["level"])
        rep = int(row["rep"])
        shot = int(row["shot"])
        shot_key = (rep, shot)
        case_level_key = (case_name, level)
        shots_per_case[case_name].add(shot_key)
        clusters_per_shot[case_level_key][shot_key] += 1
        events_by_case_level[case_level_key].append(float(row["event_count"]))
        active_by_case_level[case_level_key].append(float(row["active_detectors"]))
        diameter_by_case_level[case_level_key].append(float(row["diameter"]))
        buffer_by_case_level[case_level_key].append(float(row["buffer_bound"]))

    level_rows: list[dict[str, object]] = []
    for case_level_key in sorted(events_by_case_level.keys(), key=lambda item: (parse_distance(item[0]), item[1])):
        case_name, level = case_level_key
        events = events_by_case_level[case_level_key]
        active_sizes = active_by_case_level[case_level_key]
        diameters = diameter_by_case_level[case_level_key]
        buffers = buffer_by_case_level[case_level_key]
        shot_count = max(1, len(shots_per_case[case_name]))
        per_shot_counts = clusters_per_shot[case_level_key]
        level_rows.append(
            {
                "case": case_name,
                "distance": parse_distance(case_name),
                "level": level,
                "shots": shot_count,
                "nonempty_shots": len(per_shot_counts),
                "total_clusters": len(events),
                "avg_clusters_per_shot": f"{len(events) / shot_count:.6f}",
                "max_clusters_in_shot": max(per_shot_counts.values()) if per_shot_counts else 0,
                "total_events": f"{sum(events):.6f}",
                "mean_events_per_cluster": f"{mean(events):.6f}",
                "median_events_per_cluster": f"{median(events):.6f}",
                "p90_events_per_cluster": f"{percentile(events, 0.9):.6f}",
                "min_events_per_cluster": f"{min(events):.6f}",
                "max_events_per_cluster": f"{max(events):.6f}",
                "mean_active_detectors_per_cluster": f"{mean(active_sizes):.6f}",
                "mean_diameter_per_cluster": f"{mean(diameters):.6f}",
                "mean_buffer_bound_per_cluster": f"{mean(buffers):.6f}",
            }
        )
    return level_rows



def summarize_level_count_rows(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Summarize per-shot cluster-count distributions from level_params.csv.

    The input has one row per nonempty level in each shot. Missing levels are
    treated as cluster_count=0 so that percentiles estimate required parallel
    workers per level rather than conditional-on-nonempty counts.
    """
    if not rows:
        return [], []

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
        cm = (case, method)
        shot_key = (rep, shot)
        key = (case, method, level, rep, shot)
        shot_keys_by_case_method[cm].add(shot_key)
        levels_by_case_method[cm].add(level)
        counts[key] = int(float(row.get("cluster_count", "0") or 0))
        assigned[key] = int(float(row.get("assigned_active_detector_count", "0") or 0))
        max_active[key] = int(float(row.get("max_cluster_active_detectors", "0") or 0))
        max_events[key] = int(float(row.get("max_cluster_event_count", "0") or 0))
        max_critical[key] = int(float(row.get("max_cluster_critical_path_event_count", "0") or 0))

    summary_rows: list[dict[str, object]] = []
    histogram_rows: list[dict[str, object]] = []
    for (case, method) in sorted(shot_keys_by_case_method.keys(), key=lambda x: (parse_distance(x[0]), x[1])):
        shot_keys = sorted(shot_keys_by_case_method[(case, method)])
        max_level = max(levels_by_case_method[(case, method)]) if levels_by_case_method[(case, method)] else 0
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
            shots = len(shot_keys)
            if shots == 0:
                continue
            summary_rows.append({
                "case": case,
                "distance": parse_distance(case),
                "method": method,
                "level": level,
                "shots": shots,
                "nonempty_shots": int(sum(1 for v in cluster_values if v > 0)),
                "nonempty_fraction": f"{sum(1 for v in cluster_values if v > 0) / shots:.9f}",
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
                n = hist[cluster_count]
                histogram_rows.append({
                    "case": case,
                    "distance": parse_distance(case),
                    "method": method,
                    "level": level,
                    "cluster_count": cluster_count,
                    "shots": n,
                    "fraction": f"{n / shots:.9f}",
                })
    return summary_rows, histogram_rows


def fit_power_law(summary_rows: list[dict[str, object]], metric: str, fit_label: str, min_distance: int | None = None) -> dict[str, object] | None:
    points: list[tuple[float, float]] = []
    for row in summary_rows:
        distance = float(row["distance"])
        value = float(row[metric])
        if min_distance is not None and distance < min_distance:
            continue
        if distance > 0 and value > 0:
            points.append((distance, value))
    if len(points) < 2:
        return None
    log_x = [math.log(x) for x, _ in points]
    log_y = [math.log(y) for _, y in points]
    mean_x = sum(log_x) / len(log_x)
    mean_y = sum(log_y) / len(log_y)
    denom = sum((x - mean_x) ** 2 for x in log_x)
    if denom == 0:
        return None
    exponent = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_x, log_y)) / denom
    log_coefficient = mean_y - exponent * mean_x
    predicted = [log_coefficient + exponent * x for x in log_x]
    ss_res = sum((y - p) ** 2 for y, p in zip(log_y, predicted))
    ss_tot = sum((y - mean_y) ** 2 for y in log_y)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return {
        "metric": metric,
        "fit_range": fit_label,
        "distance_min": int(min(x for x, _ in points)),
        "distance_max": int(max(x for x, _ in points)),
        "coefficient_A": math.exp(log_coefficient),
        "exponent_b": exponent,
        "r2_log_space": r2,
        "model": f"{metric} = A * d^b",
    }


def compute_fit_rows(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    metrics = [
        "global_events_per_shot",
        "global_events_per_shot_per_round",
        "parallel_ideal_events_per_shot",
        "parallel_ideal_events_per_shot_per_round",
    ]
    fit_rows: list[dict[str, object]] = []
    for metric in metrics:
        all_fit = fit_power_law(summary_rows, metric, "all_distances")
        if all_fit is not None:
            fit_rows.append(all_fit)
        distances = sorted(int(row["distance"]) for row in summary_rows)
        if len(distances) >= 6:
            tail_min = distances[max(0, len(distances) // 3)]
            tail_fit = fit_power_law(summary_rows, metric, f"tail_d_ge_{tail_min}", tail_min)
            if tail_fit is not None:
                fit_rows.append(tail_fit)
    return fit_rows


def make_svg_plot(summary_rows: list[dict[str, object]], plot_path: Path, plot_title: str) -> None:
    distances = [int(row["distance"]) for row in summary_rows]
    global_events = [float(row["global_events_per_shot_per_round"]) for row in summary_rows]
    global_errors = [float(row.get("global_events_per_shot_per_round_std_error", 0.0)) for row in summary_rows]
    parallel_events = [float(row["parallel_ideal_events_per_shot_per_round"]) for row in summary_rows]
    parallel_errors = [float(row.get("parallel_ideal_events_per_shot_per_round_std_error", 0.0)) for row in summary_rows]
    all_y = []
    for value, err in list(zip(global_events, global_errors)) + list(zip(parallel_events, parallel_errors)):
        if value > 0:
            all_y.append(value)
            if err > 0:
                all_y.append(value + err)
                if value - err > 0:
                    all_y.append(value - err)
    if not distances or not all_y:
        raise SystemExit("cannot draw plot without positive x/y values")

    width = 940
    height = 590
    margin_left = 100
    margin_right = 35
    margin_top = 70
    margin_bottom = 85
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom

    min_log_x = math.log10(min(distances))
    max_log_x = math.log10(max(distances))
    min_log_y = math.floor(math.log10(min(all_y)))
    max_log_y = math.ceil(math.log10(max(all_y)))

    def scale_x(value: float) -> float:
        if max_log_x == min_log_x:
            return margin_left + plot_width / 2
        return margin_left + (math.log10(value) - min_log_x) / (max_log_x - min_log_x) * plot_width

    def scale_y(value: float) -> float:
        if max_log_y == min_log_y:
            return margin_top + plot_height / 2
        return margin_top + (max_log_y - math.log10(value)) / (max_log_y - min_log_y) * plot_height

    def polyline_points(xs: list[int] | list[float], ys: list[float]) -> str:
        return " ".join(f"{scale_x(x):.2f},{scale_y(y):.2f}" for x, y in zip(xs, ys) if y > 0)

    x_ticks = sorted(set(distances))
    y_ticks = [10 ** exponent for exponent in range(min_log_y, max_log_y + 1)]
    fit_rows = compute_fit_rows(summary_rows)
    global_round_fit = next((row for row in fit_rows if row["metric"] == "global_events_per_shot_per_round" and row["fit_range"] == "all_distances"), None)
    global_shot_fit = next((row for row in fit_rows if row["metric"] == "global_events_per_shot" and row["fit_range"] == "all_distances"), None)

    fit_xs: list[float] = []
    fit_ys: list[float] = []
    if global_round_fit is not None:
        a = float(global_round_fit["coefficient_A"])
        b = float(global_round_fit["exponent_b"])
        fit_xs = [min(distances) * (max(distances) / min(distances)) ** (i / 99) for i in range(100)]
        fit_ys = [a * (x ** b) for x in fit_xs]

    svg_lines: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>',
        'text { font-family: Arial, sans-serif; fill: #222; }',
        '.grid { stroke: #d9d9d9; stroke-dasharray: 4 4; }',
        '.axis { stroke: #333; stroke-width: 1.5; }',
        '.global { stroke: #0b6e4f; fill: none; stroke-width: 2.5; }',
        '.parallel { stroke: #b23a48; fill: none; stroke-width: 2.5; }',
        '.fit { stroke: #333; fill: none; stroke-width: 1.8; stroke-dasharray: 7 5; }',
        '.global-error { stroke: #0b6e4f; stroke-width: 1.2; opacity: 0.65; }',
        '.parallel-error { stroke: #b23a48; stroke-width: 1.2; opacity: 0.65; }',
        '</style>',
        f'<text x="{width / 2:.1f}" y="30" text-anchor="middle" font-size="22">{plot_title}</text>',
        f'<text x="{width / 2:.1f}" y="52" text-anchor="middle" font-size="13">Y axis is (events per shot) / rounds; error bars show SEM over shots.</text>',
        f'<line class="axis" x1="{margin_left}" y1="{margin_top + plot_height}" x2="{margin_left + plot_width}" y2="{margin_top + plot_height}" />',
        f'<line class="axis" x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{margin_top + plot_height}" />',
    ]

    for tick in y_ticks:
        y = scale_y(float(tick))
        svg_lines.append(f'<line class="grid" x1="{margin_left}" y1="{y:.2f}" x2="{margin_left + plot_width}" y2="{y:.2f}" />')
        svg_lines.append(f'<text x="{margin_left - 12}" y="{y + 5:.2f}" text-anchor="end" font-size="12">{tick:g}</text>')

    for tick in x_ticks:
        x = scale_x(float(tick))
        svg_lines.append(f'<line class="grid" x1="{x:.2f}" y1="{margin_top}" x2="{x:.2f}" y2="{margin_top + plot_height}" />')
        svg_lines.append(f'<text x="{x:.2f}" y="{margin_top + plot_height + 24}" text-anchor="middle" font-size="12">{tick}</text>')

    svg_lines.extend(
        [
            f'<text x="{margin_left + plot_width / 2:.2f}" y="{height - 22}" text-anchor="middle" font-size="15">distance d (= rounds)</text>',
            f'<text x="25" y="{margin_top + plot_height / 2:.2f}" text-anchor="middle" font-size="15" transform="rotate(-90 25 {margin_top + plot_height / 2:.2f})">events per shot per round</text>',
            f'<polyline class="global" points="{polyline_points(distances, global_events)}" />',
            f'<polyline class="parallel" points="{polyline_points(distances, parallel_events)}" />',
        ]
    )
    if fit_xs and fit_ys:
        svg_lines.append(f'<polyline class="fit" points="{polyline_points(fit_xs, fit_ys)}" />')

    def add_error_bars(xs: list[int], ys: list[float], errs: list[float], css_class: str) -> None:
        for x_value, y_value, err in zip(xs, ys, errs):
            if y_value <= 0 or err <= 0:
                continue
            y_low = max(y_value - err, min(all_y) * 0.5)
            y_high = y_value + err
            x = scale_x(float(x_value))
            y1 = scale_y(y_low)
            y2 = scale_y(y_high)
            cap = 6.0
            svg_lines.append(f'<line class="{css_class}" x1="{x:.2f}" y1="{y1:.2f}" x2="{x:.2f}" y2="{y2:.2f}" />')
            svg_lines.append(f'<line class="{css_class}" x1="{x - cap:.2f}" y1="{y1:.2f}" x2="{x + cap:.2f}" y2="{y1:.2f}" />')
            svg_lines.append(f'<line class="{css_class}" x1="{x - cap:.2f}" y1="{y2:.2f}" x2="{x + cap:.2f}" y2="{y2:.2f}" />')

    add_error_bars(distances, global_events, global_errors, "global-error")
    add_error_bars(distances, parallel_events, parallel_errors, "parallel-error")

    for x_value, y_value in zip(distances, global_events):
        if not math.isfinite(y_value) or y_value <= 0:
            continue
        svg_lines.append(
            f'<circle cx="{scale_x(x_value):.2f}" '
            f'cy="{scale_y(y_value):.2f}" '
            f'r="4.5" fill="#0b6e4f" />'
        )

    for x_value, y_value in zip(distances, parallel_events):
        if not math.isfinite(y_value) or y_value <= 0:
            continue
        x = scale_x(x_value)
        y = scale_y(y_value)
        svg_lines.append(
            f'<rect x="{x - 4.5:.2f}" y="{y - 4.5:.2f}" '
            f'width="9" height="9" fill="#b23a48" />'
        )

    legend_x = margin_left + 10
    legend_y = margin_top + 10
    svg_lines.extend(
        [
            f'<line class="global" x1="{legend_x}" y1="{legend_y}" x2="{legend_x + 34}" y2="{legend_y}" />',
            f'<circle cx="{legend_x + 17}" cy="{legend_y}" r="4.5" fill="#0b6e4f" />',
            f'<text x="{legend_x + 45}" y="{legend_y + 4}" font-size="13">global sparse blossom / round</text>',
            f'<line class="parallel" x1="{legend_x}" y1="{legend_y + 24}" x2="{legend_x + 34}" y2="{legend_y + 24}" />',
            f'<rect x="{legend_x + 12.5}" y="{legend_y + 19.5}" width="9" height="9" fill="#b23a48" />',
            f'<text x="{legend_x + 45}" y="{legend_y + 28}" font-size="13">parallel ideal level-critical / round</text>',
        ]
    )
    if global_round_fit is not None:
        exponent = float(global_round_fit["exponent_b"])
        r2 = float(global_round_fit["r2_log_space"])
        svg_lines.extend(
            [
                f'<line class="fit" x1="{legend_x}" y1="{legend_y + 48}" x2="{legend_x + 34}" y2="{legend_y + 48}" />',
                f'<text x="{legend_x + 45}" y="{legend_y + 52}" font-size="13">global / round fit: d^{exponent:.2f}, R2={r2:.4f}</text>',
            ]
        )
    if global_shot_fit is not None:
        exponent = float(global_shot_fit["exponent_b"])
        r2 = float(global_shot_fit["r2_log_space"])
        svg_lines.append(f'<text x="{legend_x}" y="{legend_y + 78}" font-size="13">raw global events/shot fit: d^{exponent:.2f}, R2={r2:.4f}</text>')
    svg_lines.append('</svg>')

    plot_path.parent.mkdir(parents=True, exist_ok=True)
    plot_path.write_text("\n".join(svg_lines), encoding="utf-8")

def main() -> None:
    args = parse_args()
    bench_rows = read_rows(args.bench_csv)
    cluster_rows = read_rows(args.cluster_csv)
    level_param_rows = read_rows(args.level_params_csv) if args.level_params_csv is not None else []
    if not bench_rows:
        raise SystemExit(f"no rows found in {args.bench_csv}")

    summary_rows = summarize_bench_rows(bench_rows)
    level_rows = summarize_cluster_rows(cluster_rows) if cluster_rows else []
    level_count_summary_rows, level_count_histogram_rows = summarize_level_count_rows(level_param_rows)
    if not summary_rows:
        raise SystemExit("no comparable global/parallel rows found in benchmark csv")

    write_csv(
        args.summary_csv,
        [
            "case",
            "distance",
            "parallel_mode",
            "global_events_per_shot",
            "global_events_std_error_per_shot",
            "parallel_ideal_events_per_shot",
            "parallel_ideal_events_std_error_per_shot",
            "parallel_total_cluster_events_per_shot",
            "parallel_total_cluster_events_std_error_per_shot",
            "global_events_per_shot_per_round",
            "global_events_per_shot_per_round_std_error",
            "parallel_ideal_events_per_shot_per_round",
            "parallel_ideal_events_per_shot_per_round_std_error",
            "parallel_total_cluster_events_per_shot_per_round",
            "parallel_total_cluster_events_per_shot_per_round_std_error",
            "parallel_nonempty_event_levels_per_shot",
            "parallel_critical_path_algorithmic_time_per_shot",
            "parallel_critical_path_algorithmic_time_std_error_per_shot",
            "avg_clusters",
            "avg_max_level",
            "ideal_event_count_speedup_global_over_parallel",
            "mistakes",
            "exceptions",
        ],
        summary_rows,
    )
    write_csv(
        args.level_stats_csv,
        [
            "case",
            "distance",
            "level",
            "shots",
            "nonempty_shots",
            "total_clusters",
            "avg_clusters_per_shot",
            "max_clusters_in_shot",
            "total_events",
            "mean_events_per_cluster",
            "median_events_per_cluster",
            "p90_events_per_cluster",
            "min_events_per_cluster",
            "max_events_per_cluster",
            "mean_active_detectors_per_cluster",
            "mean_diameter_per_cluster",
            "mean_buffer_bound_per_cluster",
        ],
        level_rows,
    )
    if args.level_count_summary_csv is not None:
        write_csv(
            args.level_count_summary_csv,
            [
                "case",
                "distance",
                "method",
                "level",
                "shots",
                "nonempty_shots",
                "nonempty_fraction",
                "mean_cluster_count",
                "median_cluster_count",
                "p90_cluster_count",
                "p95_cluster_count",
                "p99_cluster_count",
                "max_cluster_count",
                "mean_assigned_active_detector_count",
                "p95_assigned_active_detector_count",
                "max_assigned_active_detector_count",
                "mean_max_cluster_active_detectors",
                "p95_max_cluster_active_detectors",
                "max_cluster_active_detectors",
                "mean_max_cluster_event_count",
                "max_cluster_event_count",
                "mean_max_cluster_critical_path_event_count",
                "max_cluster_critical_path_event_count",
            ],
            level_count_summary_rows,
        )
    if args.level_count_histogram_csv is not None:
        write_csv(
            args.level_count_histogram_csv,
            [
                "case",
                "distance",
                "method",
                "level",
                "cluster_count",
                "shots",
                "fraction",
            ],
            level_count_histogram_rows,
        )

    write_csv(
        args.plot_data_csv,
        [
            "case",
            "distance",
            "global_events_per_shot",
            "global_events_std_error_per_shot",
            "parallel_ideal_events_per_shot",
            "parallel_ideal_events_std_error_per_shot",
            "parallel_total_cluster_events_per_shot",
            "parallel_total_cluster_events_std_error_per_shot",
            "global_events_per_shot_per_round",
            "global_events_per_shot_per_round_std_error",
            "parallel_ideal_events_per_shot_per_round",
            "parallel_ideal_events_per_shot_per_round_std_error",
            "parallel_total_cluster_events_per_shot_per_round",
            "parallel_total_cluster_events_per_shot_per_round_std_error",
            "ideal_event_count_speedup_global_over_parallel",
        ],
        summary_rows,
    )
    fit_rows = compute_fit_rows(summary_rows)
    if args.fit_csv is not None:
        write_csv(
            args.fit_csv,
            [
                "metric",
                "fit_range",
                "distance_min",
                "distance_max",
                "coefficient_A",
                "exponent_b",
                "r2_log_space",
                "model",
            ],
            [
                {
                    **row,
                    "coefficient_A": f"{float(row['coefficient_A']):.9g}",
                    "exponent_b": f"{float(row['exponent_b']):.9f}",
                    "r2_log_space": f"{float(row['r2_log_space']):.9f}",
                }
                for row in fit_rows
            ],
        )
    make_svg_plot(summary_rows, args.plot_svg, args.plot_title)


if __name__ == "__main__":
    main()
