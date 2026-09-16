#!/usr/bin/env python3
"""
Generate sparse-blossom benchmark case files from Stim surface-code circuits.

Default target:
    d=17, rounds=17, p=0.001, shots=256, seed=217

Noise policy (default):
    - 2-qubit Clifford error = p
    - 1-qubit Clifford error = p / 10
    - |0>, |+> preparation error = p
    - Z/X measurement error = p
    - idle error = 0

Stim exposes one ``after_clifford_depolarization`` parameter for both 1Q and
2Q Cliffords.  We therefore generate Clifford noise at p and recursively
retune the generated DEPOLARIZE1 instructions to p/10.  The explicit idle
pass is skipped when idle_factor=0.

The output case format matches src/local_parallel_bench.cc:
    <num_detectors> <num_observables> <num_edges> <num_shots> <rounds> <distance> <noise>
    e <u> <v|-1> <weight> <probability> <num_obs_flips> [obs...]
    s <obs_mask> <num_hits> [detector ids...]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, Sequence

import stim


def insert_idle_depolarization(circuit: stim.Circuit, p: float) -> stim.Circuit:
    """
    各TICKごとに idle qubit に DEPOLARIZE1(p) を挿入する。
    REPEAT ブロックも再帰的に処理する。
    MR の後に X_ERROR がある場合は、その qubit に idle depolarization を2回追加。
    """

    def process_block(block: stim.Circuit, all_qubits: set[int], repeat_flag: bool) -> stim.Circuit:
        flag = repeat_flag
        new_block = stim.Circuit()
        tick_qubits = set()
        last_was_mr = False
        last_mr_qubits = set()

        for inst in block:
            name = inst.name

            # QUBIT_COORDS はそのままコピー
            if name == "QUBIT_COORDS":
                new_block.append(inst)
                continue

            # REPEATブロックは再帰的に処理
            if name == "REPEAT":
                repeat_count = inst.repeat_count
                processed_body = process_block(inst.body_copy(), all_qubits, 1)
                new_block.append(stim.CircuitRepeatBlock(repeat_count, processed_body))
                continue

            # TICK の前に idle depolarization を入れる
            if name == "TICK":
                if flag == 1:
                    flag = 0
                else:
                    idle = sorted(all_qubits - tick_qubits)
                    if idle:
                        new_block.append(stim.CircuitInstruction("DEPOLARIZE1", idle, [p]))
                    new_block.append(inst)
                    tick_qubits.clear()
                    last_was_mr = False
                    last_mr_qubits.clear()
                    continue

            # MR → 次の X_ERROR の特別処理のため記録
            if name == "MR":
                mr_qubits = set()
                for g in inst.target_groups():
                    for t in g:
                        mr_qubits.add(t.value)
                        tick_qubits.add(t.value)
                new_block.append(inst)
                last_was_mr = True
                last_mr_qubits = mr_qubits
                continue

            # MRの直後のX_ERRORの特別処理
            if last_was_mr and name == "X_ERROR":
                err_qubits = {t.value for g in inst.target_groups() for t in g}
                if err_qubits == last_mr_qubits:
                    # MR直後のX_ERRORなので2回idle depol追加（後でMR idleに対応）
                    new_block.append(inst)
                    new_block.append(stim.CircuitInstruction("DEPOLARIZE1", sorted(all_qubits - err_qubits), [p]))
                    new_block.append(stim.CircuitInstruction("DEPOLARIZE1", sorted(all_qubits - err_qubits), [p]))
                    last_was_mr = False
                    continue

            # 通常命令
            new_block.append(inst)
            for g in inst.target_groups():
                for t in g:
                    tick_qubits.add(t.value)
            last_was_mr = False
            last_mr_qubits.clear()

        # ブロック末尾に idle depolarization を追加（最後のTICK後）
        idle = sorted(all_qubits - tick_qubits)
        if idle:
            new_block.append(stim.CircuitInstruction("DEPOLARIZE1", idle, [p]))

        return new_block

    # 全 qubit index を抽出
    all_qubits = set()
    for inst in circuit:
        if inst.name == "QUBIT_COORDS":
            all_qubits.add(inst.target_groups()[0][0].value)

    # 処理本体
    return process_block(circuit, all_qubits, 0)


def retune_single_qubit_clifford_noise(circuit: stim.Circuit, p1q: float) -> stim.Circuit:
    """Replace generated 1Q Clifford DEPOLARIZE1 probabilities, including inside REPEAT blocks."""
    out = stim.Circuit()
    for inst in circuit:
        if inst.name == "REPEAT":
            body = retune_single_qubit_clifford_noise(inst.body_copy(), p1q)
            out.append(stim.CircuitRepeatBlock(inst.repeat_count, body))
        elif inst.name == "DEPOLARIZE1":
            targets = [t for group in inst.target_groups() for t in group]
            out.append(stim.CircuitInstruction("DEPOLARIZE1", targets, [p1q]))
        else:
            out.append(inst)
    return out


def make_base_surface_code_circuit(
    *,
    distance: int,
    rounds: int,
    p: float,
    p1q_factor: float,
    before_round_data_factor: float,
    task: str,
) -> stim.Circuit:
    """Create the target surface-code circuit before optional explicit idle noise."""
    circuit = stim.Circuit.generated(
        task,
        distance=distance,
        rounds=rounds,
        # Stim uses this one parameter for both 1Q and 2Q Cliffords.
        # DEPOLARIZE1 is retuned below, while DEPOLARIZE2 stays at p.
        after_clifford_depolarization=p,
        after_reset_flip_probability=p,
        before_measure_flip_probability=p,
        before_round_data_depolarization=p * before_round_data_factor,
    )
    return retune_single_qubit_clifford_noise(circuit, p * p1q_factor)


def make_surface_code_circuit(
    *,
    distance: int,
    rounds: int,
    p: float,
    p1q_factor: float,
    idle_factor: float,
    before_round_data_factor: float,
    task: str,
) -> stim.Circuit:
    """Create a Stim circuit with 2Q=p, 1Q=p*p1q_factor, prep/meas=p, idle=p*idle_factor."""
    base_circuit = make_base_surface_code_circuit(
        distance=distance,
        rounds=rounds,
        p=p,
        p1q_factor=p1q_factor,
        before_round_data_factor=before_round_data_factor,
        task=task,
    )
    p_idle = p * idle_factor
    if p_idle == 0:
        return base_circuit
    return insert_idle_depolarization(base_circuit, p_idle)


def dem_target_value(target) -> int:
    """Compatibility wrapper for Stim DemTarget integer values."""
    if hasattr(target, "val"):
        return int(target.val)
    return int(target.value)


def iter_graphlike_dem_errors(dem: stim.DetectorErrorModel):
    """
    Yield graphlike DEM errors as (u, v, probability, observable_ids).

    v == -1 means a boundary edge.  Errors with separators or more than two
    detector targets are rejected because local_parallel_bench.cc expects a
    graphlike matching graph.
    """
    for inst in dem:
        inst_type = getattr(inst, "type", getattr(inst, "name", None))
        if inst_type == "detector" or inst_type == "logical_observable":
            continue
        if inst_type == "repeat":
            raise ValueError("DEM still contains a repeat block; call detector_error_model(flatten_loops=True).")
        if inst_type != "error":
            continue

        args = list(inst.args_copy())
        if len(args) != 1:
            raise ValueError(f"Unsupported DEM error arguments: {args!r}")
        probability = float(args[0])
        # Stim uses separator targets (^) to mark a decomposed error made of
        # graphlike pieces.  The local benchmark consumes a simple matching graph
        # without correlation metadata, so split such an instruction into the
        # graphlike pieces and assign the DEM probability to each piece.
        pieces: list[tuple[list[int], list[int]]] = [([], [])]
        for target in inst.targets_copy():
            if target.is_separator():
                if pieces[-1][0] or pieces[-1][1]:
                    pieces.append(([], []))
                continue
            if target.is_relative_detector_id():
                pieces[-1][0].append(dem_target_value(target))
            elif target.is_logical_observable_id():
                pieces[-1][1].append(dem_target_value(target))
            else:
                raise ValueError(f"Unsupported DEM target: {target!r}")

        for detectors, observables in pieces:
            observables = sorted(set(observables))
            if len(detectors) == 1:
                yield detectors[0], -1, probability, observables
            elif len(detectors) == 2:
                u, v = sorted(detectors)
                yield u, v, probability, observables
            elif len(detectors) == 0 and not observables:
                # Purely silent probability term. It is irrelevant to matching.
                continue
            else:
                raise ValueError(
                    f"Non-graphlike DEM piece with {len(detectors)} detectors and observables {observables}."
                )


def matching_weight_from_probability(p: float) -> float:
    """Return log((1-p)/p), matching the weight convention used by UserGraph."""
    if not 0 < p < 1:
        if p == 0:
            return math.inf
        if p == 1:
            return -math.inf
        raise ValueError(f"Probability must be in [0, 1], got {p!r}")
    return math.log((1 - p) / p)


def obs_bits_to_mask(bits: Sequence[bool]) -> int:
    if len(bits) > 64:
        raise ValueError(f"local_parallel_bench.cc uses uint64_t obs_mask, but this circuit has {len(bits)} observables.")
    mask = 0
    for k, bit in enumerate(bits):
        if bit:
            mask |= 1 << k
    return mask


def write_case_file(
    *,
    circuit: stim.Circuit,
    output_path: Path,
    distance: int,
    rounds: int,
    p: float,
    shots: int,
    seed: int,
) -> None:
    """Write a benchmark case file consumed by src/local_parallel_bench.cc."""
    dem = circuit.detector_error_model(decompose_errors=True, flatten_loops=True)
    edges = list(iter_graphlike_dem_errors(dem))

    sampler = circuit.compile_detector_sampler(seed=seed)
    detection_events, observable_flips = sampler.sample(shots=shots, separate_observables=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        f.write(
            f"{dem.num_detectors} {dem.num_observables} {len(edges)} "
            f"{shots} {rounds} {distance} {p}\n"
        )
        for u, v, probability, observables in edges:
            weight = matching_weight_from_probability(probability)
            obs_suffix = "" if not observables else " " + " ".join(str(o) for o in observables)
            f.write(f"e {u} {v} {weight:.17g} {probability:.17g} {len(observables)}{obs_suffix}\n")

        for shot_index in range(shots):
            hits = [str(k) for k, bit in enumerate(detection_events[shot_index]) if bit]
            obs_mask = obs_bits_to_mask(observable_flips[shot_index])
            hit_suffix = "" if not hits else " " + " ".join(hits)
            f.write(f"s {obs_mask} {len(hits)}{hit_suffix}\n")


def default_output_path(distance: int, rounds: int, p: float, shots: int, seed: int) -> Path:
    p_text = f"{p:.4f}"
    return Path("cases") / f"case_r{rounds}_d{distance}_p{p_text}_s{shots}_seed{seed}_full_idle_with_brd.txt"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance", type=int, default=17)
    parser.add_argument("--rounds", type=int, default=None, help="Defaults to --distance.")
    parser.add_argument("--p", type=float, default=0.001, help="Base error rate: 2Q, preparation, and measurement.")
    parser.add_argument("--p1q-factor", type=float, default=0.1, help="1Q Clifford error rate as a factor of --p (default: 0.1).")
    parser.add_argument("--idle-factor", type=float, default=0.0, help="Explicit idle error rate as a factor of --p (default: 0).")
    parser.add_argument("--before-round-data-factor", type=float, default=0.0, help="Stim before_round_data_depolarization as a factor of --p (uniform compatibility: 1.0).")
    parser.add_argument("--shots", type=int, default=256)
    parser.add_argument("--seed", type=int, default=217)
    parser.add_argument("--task", default="surface_code:rotated_memory_x")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--write-stim", type=Path, default=None, help="Optional path for the generated .stim circuit.")
    parser.add_argument("--write-dem", type=Path, default=None, help="Optional path for the flattened/decomposed .dem model.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rounds = args.distance if args.rounds is None else args.rounds
    output_path = args.out or default_output_path(args.distance, rounds, args.p, args.shots, args.seed)

    circuit = make_surface_code_circuit(
        distance=args.distance,
        rounds=rounds,
        p=args.p,
        p1q_factor=args.p1q_factor,
        idle_factor=args.idle_factor,
        before_round_data_factor=args.before_round_data_factor,
        task=args.task,
    )

    if args.write_stim is not None:
        args.write_stim.parent.mkdir(parents=True, exist_ok=True)
        args.write_stim.write_text(str(circuit), encoding="utf-8")

    if args.write_dem is not None:
        dem = circuit.detector_error_model(decompose_errors=True, flatten_loops=True)
        args.write_dem.parent.mkdir(parents=True, exist_ok=True)
        args.write_dem.write_text(str(dem), encoding="utf-8")

    write_case_file(
        circuit=circuit,
        output_path=output_path,
        distance=args.distance,
        rounds=rounds,
        p=args.p,
        shots=args.shots,
        seed=args.seed,
    )
    print(output_path)


if __name__ == "__main__":
    main()
