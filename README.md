# Parallel Sparse Blossom: event-count reproduction code

This repository contains the code and numerical data used for the event-count experiments in the accompanying manuscript on parallel sparse-blossom decoding.

The public package contains the **fixed parallel sparse-blossom** experiment used in the current manuscript, together with the non-parallel global sparse-blossom baseline used for comparison.

## Repository layout

```text
parallel_sparse_blossom/
├── README.md
├── LICENSE-PyMatching
├── requirements.txt
├── data/
│   ├── README.md
│   └── figure4_event_counts.csv     # numerical values underlying Fig. 4
├── src/
│   ├── local_parallel_bench.cc
│   ├── stim.h
│   └── pymatching/                  # modified sparse-blossom / PyMatching sources
├── scripts/
│   ├── build.sh                     # compile the C++ benchmark
│   ├── gen_cases.sh                 # generate Stim benchmark cases
│   ├── run_compare.sh               # main fixed-vs-global experiment
│   ├── run_sweep.sh                 # internal distance-sweep runner
│   ├── analyze.py                   # summarize one sweep
│   ├── plot_compare.py              # global-vs-fixed comparison plot
│   ├── merge_level_counts.py        # level/core-count summaries
│   └── lib.sh
└── tools/
    └── generate_stim_case_with_idle.py
```

Generated directories (`build/`, `bin/`, `cases/`, and `results/`) are intentionally not included in the archive.

The compact `data/figure4_event_counts.csv` file contains the numerical values underlying Fig. 4 of the manuscript.

## Requirements

- Linux or another Unix-like environment with Bash
- Python 3
- A C++20 compiler (`g++` by default)
- Python package `stim`

Install the Python dependency with

```bash
python3 -m pip install -r requirements.txt
```

The build script uses `g++ -std=c++20 -O3 -pthread` by default. You can override the compiler or optimization flag with `CXX` and `OPT`.

## Noise models used in the numerical experiments

The case generator supports both circuit-level noise models used in Fig. 4 through `--noise-model`:

- `physical` (default): two-qubit Clifford error `p`, one-qubit Clifford error `p/10`, preparation and measurement error `p`, and idle error `0`;
- `uniform`: two-qubit Clifford, one-qubit Clifford, preparation, measurement, and idle error rates are all parameterized by `p` as described in the manuscript.

For the generated depolarizing channels, the total one-qubit error probability is the value above (distributed uniformly over `X`, `Y`, and `Z`), and the total two-qubit error probability is `p` (distributed uniformly over the 15 nonidentity two-qubit Paulis).

The historical filename suffix `_full_idle_with_brd.txt` is retained for compatibility with existing case files.

## Build

From the repository root:

```bash
bash ./scripts/build.sh
```

This creates

```text
bin/parallel_eventcount_bench
```

## Generate benchmark cases

For example, to generate the `p = 10^-5`, `d = 9, 13, ..., 49` cases with 256 shots:

```bash
bash ./scripts/gen_cases.sh \
  --case-dir cases/default_p1e5_d9_d49_s256 \
  --distances '9:49:4' \
  --p 0.00001 \
  --noise-model physical \
  --shots 256 \
  --seed-base 200
```

The distance syntax `9:49:4` means start at 9, stop at 49, and increment by 4.

## Run the main experiment

### Complete reproduction command

`run_compare.sh` can generate the cases and build the benchmark automatically. A normal foreground command corresponding to the paper-style `p = 10^-5` run is:

```bash
bash ./scripts/run_compare.sh \
  --case-dir cases/default_p1e5_d9_d49_s256 \
  --case-template 'case_r{distance}_d{distance}_p{p4}_s256_seed{seed}_full_idle_with_brd.txt' \
  --p 0.00001 \
  --noise-model physical \
  --workers 8 \
  --distances '9:49:4' \
  --max-shots 256 \
  --reps 1 \
  --parameter-schedule \
  --phi-sequence-floor 0.01 \
  --phi-sequence-q 0.1
```

The current paper parameter schedule is the default, so `--parameter-schedule` is included above mainly for clarity.

To generate the corresponding uniform-noise data, replace `--noise-model physical` by `--noise-model uniform`. The same choice can be made for each of the three physical error rates `1e-3`, `1e-4`, and `1e-5`.

### Command used when cases and the binary already exist

The long runs used pre-generated cases and a pre-built binary. The corresponding foreground command is:

```bash
bash ./scripts/run_compare.sh \
  --case-dir cases/default_p1e5_d9_d49_s256 \
  --case-template 'case_r{distance}_d{distance}_p0.0000_s256_seed{seed}_full_idle_with_brd.txt' \
  --p 0.00001 \
  --workers 8 \
  --distances '9:49:4' \
  --max-shots 256 \
  --reps 1 \
  --parameter-schedule \
  --phi-sequence-floor 0.01 \
  --phi-sequence-q 0.1 \
  --no-generate-cases \
  --no-build
```

This is the same command structure as the background run used for the experiment, but without `nohup`, shell redirection, or `&`.

The explicit `--p 0.00001` keeps the recorded run metadata consistent with the pre-generated cases. The historical filename uses four-decimal formatting, so `p = 10^-5` appears as `p0.0000`.

## What `run_compare.sh` measures

For each selected distance, the script runs the same case through:

1. the global sparse-blossom baseline; and
2. the fixed parallel sparse-blossom event-count implementation.

The fixed hierarchy uses the paper's parameter schedule by default. The benchmark also checks the observable prediction against the global sparse-blossom result and reports `mistakes` and `exceptions` in the output CSVs.

## Main output files

Unless `--out-dir` is specified, a timestamped directory is created under `results/`:

```text
results/compare_<timestamp>/
├── compare.csv
├── compare.svg
├── config.env
├── level_cluster_count_summary.csv
├── level_cluster_count_histogram.csv
└── fixed/
    ├── bench_rows.csv
    ├── bench_stderr.log
    ├── eventcount_summary.csv
    ├── level_params.csv
    └── run_config.env
```

Additional debug CSVs are produced only when the corresponding logging options are enabled.

Useful options include:

```text
--assignment-log 1
--cluster-log 1
--mistake-log 1
--sweep-plot 1
--fixed-bounds D1:B1,D2:B2,...
```

Run

```bash
bash ./scripts/run_compare.sh --help
```

for the complete option list.

## Running the lower-level sweep directly

For most reproduction work, use `run_compare.sh`. `scripts/run_sweep.sh` is the lower-level runner called by it and is useful when directly controlling generated cases, environment variables, or debug output.

## Numerical data underlying Figure 4

The repository includes `data/figure4_event_counts.csv`, a compact table containing all points plotted in Fig. 4: three physical error rates (`1e-3`, `1e-4`, `1e-5`) for both the uniform and physically motivated noise models, at distances `d = 9, 13, ..., 49`. Each point is based on 256 shots.

The table contains the global and parallel event counts per shot per syndrome-extraction round together with their standard errors. The dotted `d^2` curves shown in the figure are analytic visual guides and are not stored as data. See `data/README.md` for the column definitions.

## Implementation notes

The C++ implementation contains a modified sparse-blossom/PyMatching code path for hierarchical processing clusters, state import/export, event-count accounting, and final observable aggregation. The public fixed-only version includes the observable-extraction handling for touched inherited active-detector sources and the paper Algorithm 3 parameter schedule, including the `2 w_max` additive term in the next-level diameter bound.

The package deliberately omits temporary parameter-launch scripts, background-job helpers, patches, logs, generated cases, and bulky intermediate/debug result files. The compact numerical table underlying Fig. 4 is retained under `data/`.

## Third-party code

Files under `src/pymatching/` are derived from [**PyMatching**](https://github.com/oscarhiggott/PyMatching) presented by Oscar Higgott, Craig Gidney and other contributors, and retain their original copyright and Apache-2.0 license headers. A copy of the PyMatching Apache-2.0 license is included as `LICENSE-PyMatching`. No changes have been made to the files copied from PyMatching.

Stim is used as a Python dependency for case generation and is not vendored in this repository.
