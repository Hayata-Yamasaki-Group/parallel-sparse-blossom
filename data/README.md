# Numerical data used in Figure 4

`figure4_event_counts.csv` contains the numerical values underlying Figure 4 of the manuscript.

The file combines the six simulation conditions shown in the figure:

- uniform circuit-level i.i.d. depolarizing noise at `p = 1e-3`, `1e-4`, and `1e-5`;
- physically motivated circuit-level noise at `p = 1e-3`, `1e-4`, and `1e-5`.

For every condition, the code distances are `d = 9, 13, ..., 49`, and each point is based on 256 shots.

## Columns

- `noise_model`: `uniform` or `physically_motivated`.
- `physical_error_rate`: nominal physical error rate `p`.
- `distance`: surface-code distance `d`; the number of syndrome-extraction rounds is also `d`.
- `shots`: number of simulated shots.
- `global_events_per_shot_per_round`: mean event count of global sparse blossom, divided by the number of rounds.
- `global_events_per_shot_per_round_std_error`: standard error of that mean.
- `parallel_events_per_shot_per_round`: mean critical-path event count of the fixed parallel sparse-blossom implementation, divided by the number of rounds.
- `parallel_events_per_shot_per_round_std_error`: standard error of that mean.

The dotted `d^2` curves in Figure 4 are visual guides and therefore are not stored as numerical data.

The CSV was assembled directly from the `compare.csv` files produced by the numerical runs; no plotted numerical values were refitted or modified during this cleanup.
