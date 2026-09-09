#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

CASE_DISTANCE_RE = re.compile(r"_d(?P<distance>\d+)_")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot global and fixed-clustering event counts.")
    p.add_argument("--fixed-summary-csv", required=True, type=Path)
    p.add_argument("--comparison-csv", required=True, type=Path)
    p.add_argument("--plot-svg", required=True, type=Path)
    p.add_argument("--fixed-label", default="fixed ideal")
    p.add_argument("--plot-title", default="Global vs Fixed Event Counts")
    return p.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def parse_distance(case_name: str) -> int:
    m = CASE_DISTANCE_RE.search(case_name)
    if not m:
        raise ValueError(f"could not parse distance from case name: {case_name}")
    return int(m.group("distance"))


def f(row: dict[str, str], key: str) -> float:
    raw = row.get(key, "")
    return float(raw) if raw not in ("", None) else 0.0


def merge_rows(fixed_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for x in sorted(fixed_rows, key=lambda r: int(r.get("distance") or parse_distance(r.get("case", "")))):
        d = int(x.get("distance") or parse_distance(x.get("case", "")))
        global_per_shot = f(x, "global_events_per_shot")
        fixed_per_shot = f(x, "parallel_ideal_events_per_shot")
        rows.append({
            "distance": d,
            "global_case": x.get("case", ""),
            "fixed_case": x.get("case", ""),
            "global_events_per_shot": f"{global_per_shot:.9f}",
            "global_events_std_error_per_shot": f"{f(x, 'global_events_std_error_per_shot'):.9f}",
            "global_events_per_shot_per_round": f"{f(x, 'global_events_per_shot_per_round'):.12f}",
            "global_events_per_shot_per_round_std_error": f"{f(x, 'global_events_per_shot_per_round_std_error'):.12f}",
            "fixed_ideal_events_per_shot": f"{fixed_per_shot:.9f}",
            "fixed_ideal_events_std_error_per_shot": f"{f(x, 'parallel_ideal_events_std_error_per_shot'):.9f}",
            "fixed_ideal_events_per_shot_per_round": f"{f(x, 'parallel_ideal_events_per_shot_per_round'):.12f}",
            "fixed_ideal_events_per_shot_per_round_std_error": f"{f(x, 'parallel_ideal_events_per_shot_per_round_std_error'):.12f}",
            "speedup_global_over_fixed": f"{(global_per_shot / fixed_per_shot) if fixed_per_shot > 0 else 0.0:.9f}",
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise SystemExit("no rows in fixed summary CSV file")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fobj:
        writer = csv.DictWriter(fobj, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def fit_power_law(points: list[tuple[float, float]]) -> tuple[float, float, float] | None:
    if len(points) < 2:
        return None
    lx = [math.log(x) for x, _ in points]
    ly = [math.log(y) for _, y in points]
    mx = sum(lx) / len(lx)
    my = sum(ly) / len(ly)
    denom = sum((x - mx) ** 2 for x in lx)
    if denom == 0:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(lx, ly)) / denom
    a_log = my - b * mx
    pred = [a_log + b * x for x in lx]
    ss_res = sum((y - p) ** 2 for y, p in zip(ly, pred))
    ss_tot = sum((y - my) ** 2 for y in ly)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return math.exp(a_log), b, r2


def make_svg(rows: list[dict[str, object]], out: Path, title: str, fixed_label: str) -> None:
    d = [int(r["distance"]) for r in rows]
    gy = [float(r["global_events_per_shot_per_round"]) for r in rows]
    ge = [float(r["global_events_per_shot_per_round_std_error"]) for r in rows]
    fy = [float(r["fixed_ideal_events_per_shot_per_round"]) for r in rows]
    fe = [float(r["fixed_ideal_events_per_shot_per_round_std_error"]) for r in rows]
    all_y = [v for v in gy + fy if v > 0]
    all_y += [v + e for vals, errs in [(gy, ge), (fy, fe)] for v, e in zip(vals, errs) if v > 0 and e > 0]
    if not d or not all_y:
        raise SystemExit("cannot draw comparison plot without positive values")

    width, height = 980, 620
    ml, mr, mt, mb = 105, 35, 70, 90
    pw, ph = width - ml - mr, height - mt - mb
    min_log_x, max_log_x = math.log10(min(d)), math.log10(max(d))
    min_log_y = math.floor(math.log10(min(all_y)))
    max_log_y = math.ceil(math.log10(max(all_y)))

    def sx(v: float) -> float:
        return ml + (math.log10(v) - min_log_x) / (max_log_x - min_log_x) * pw if max_log_x != min_log_x else ml + pw / 2

    def sy(v: float) -> float:
        return mt + (max_log_y - math.log10(v)) / (max_log_y - min_log_y) * ph if max_log_y != min_log_y else mt + ph / 2

    def poly(xs: list[int], ys: list[float]) -> str:
        return " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(xs, ys) if y > 0)

    def esc(s: str) -> str:
        return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text { font-family: Arial, sans-serif; fill: #222; }",
        ".grid { stroke: #d9d9d9; stroke-dasharray: 4 4; }",
        ".axis { stroke: #333; stroke-width: 1.5; }",
        ".global { stroke: #0b6e4f; fill: none; stroke-width: 2.5; }",
        ".fixed { stroke: #2b6cb0; fill: none; stroke-width: 2.5; }",
        ".errg { stroke: #0b6e4f; stroke-width: 1.1; opacity: 0.65; }",
        ".errf { stroke: #2b6cb0; stroke-width: 1.1; opacity: 0.65; }",
        "</style>",
        f'<text x="{width/2:.1f}" y="30" text-anchor="middle" font-size="22">{esc(title)}</text>',
        f'<text x="{width/2:.1f}" y="52" text-anchor="middle" font-size="13">Y axis is (events per shot) / rounds; error bars show SEM over shots.</text>',
        f'<line class="axis" x1="{ml}" y1="{mt + ph}" x2="{ml + pw}" y2="{mt + ph}" />',
        f'<line class="axis" x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt + ph}" />',
    ]
    for exp in range(min_log_y, max_log_y + 1):
        tick = 10 ** exp
        y = sy(float(tick))
        svg += [f'<line class="grid" x1="{ml}" y1="{y:.2f}" x2="{ml + pw}" y2="{y:.2f}" />',
                f'<text x="{ml - 12}" y="{y + 5:.2f}" text-anchor="end" font-size="12">{tick:g}</text>']
    for tick in sorted(set(d)):
        x = sx(float(tick))
        svg += [f'<line class="grid" x1="{x:.2f}" y1="{mt}" x2="{x:.2f}" y2="{mt + ph}" />',
                f'<text x="{x:.2f}" y="{mt + ph + 24}" text-anchor="middle" font-size="12">{tick}</text>']
    svg += [
        f'<text x="{ml + pw/2:.2f}" y="{height - 22}" text-anchor="middle" font-size="15">distance d (= rounds)</text>',
        f'<text x="25" y="{mt + ph/2:.2f}" text-anchor="middle" font-size="15" transform="rotate(-90 25 {mt + ph/2:.2f})">events per shot per round</text>',
        f'<polyline class="global" points="{poly(d, gy)}" />',
        f'<polyline class="fixed" points="{poly(d, fy)}" />',
    ]

    def err(xs: list[int], ys: list[float], es: list[float], cls: str) -> None:
        cap = 6.0
        floor = min(all_y) * 0.5
        for xval, yval, e in zip(xs, ys, es):
            if yval <= 0 or e <= 0:
                continue
            x = sx(float(xval))
            y1 = sy(max(yval - e, floor))
            y2 = sy(yval + e)
            svg.append(f'<line class="{cls}" x1="{x:.2f}" y1="{y1:.2f}" x2="{x:.2f}" y2="{y2:.2f}" />')
            svg.append(f'<line class="{cls}" x1="{x-cap:.2f}" y1="{y1:.2f}" x2="{x+cap:.2f}" y2="{y1:.2f}" />')
            svg.append(f'<line class="{cls}" x1="{x-cap:.2f}" y1="{y2:.2f}" x2="{x+cap:.2f}" y2="{y2:.2f}" />')
    err(d, gy, ge, "errg")
    err(d, fy, fe, "errf")
    for x, y in zip(d, gy):
        if y > 0:
            svg.append(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="4.4" fill="#0b6e4f" />')
    for x, y in zip(d, fy):
        if y > 0:
            svg.append(f'<rect x="{sx(x)-4.4:.2f}" y="{sy(y)-4.4:.2f}" width="8.8" height="8.8" fill="#2b6cb0" />')

    lx, ly = ml + 10, mt + 10
    svg += [
        f'<line class="global" x1="{lx}" y1="{ly}" x2="{lx+34}" y2="{ly}" /><circle cx="{lx+17}" cy="{ly}" r="4.4" fill="#0b6e4f" />',
        f'<text x="{lx+45}" y="{ly+4}" font-size="13">global sparse blossom / round</text>',
        f'<line class="fixed" x1="{lx}" y1="{ly+24}" x2="{lx+34}" y2="{ly+24}" /><rect x="{lx+12.6}" y="{ly+19.6}" width="8.8" height="8.8" fill="#2b6cb0" />',
        f'<text x="{lx+45}" y="{ly+28}" font-size="13">{esc(fixed_label)} / round</text>',
    ]
    ytxt = ly + 60
    for label, ys in [("global", gy), (fixed_label, fy)]:
        fit = fit_power_law([(float(x), float(y)) for x, y in zip(d, ys) if y > 0])
        if fit:
            _a, b, r2 = fit
            svg.append(f'<text x="{lx}" y="{ytxt}" font-size="12">fit {esc(label)}: A*d^{b:.2f}, R2={r2:.4f}</text>')
            ytxt += 18
    svg.append("</svg>")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(svg), encoding="utf-8")


def main() -> None:
    args = parse_args()
    rows = merge_rows(read_rows(args.fixed_summary_csv))
    write_csv(args.comparison_csv, rows)
    make_svg(rows, args.plot_svg, args.plot_title, args.fixed_label)


if __name__ == "__main__":
    main()
