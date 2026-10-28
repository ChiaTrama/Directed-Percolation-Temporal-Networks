#!/usr/bin/env python3
"""Small-scale reproduction of subcritical cluster-mass growth.

The model follows Badie-Modiri et al., PRE 105, 054313 (2022): a random
k-regular base graph, independent Poisson activations on every edge, and a
reduced temporal event graph.  It is designed as a transparent laptop-scale
illustration, not as a replacement for the authors' high-statistics code.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-complex-trovato")

import matplotlib.pyplot as plt
import numpy as np


ANALYTICAL_THRESHOLD_K9 = 0.08559
EMPIRICAL_THRESHOLD_K9 = 0.08808


@dataclass(frozen=True)
class Config:
    nodes: int
    degree: int
    duration: float
    mean_dt: float
    realisations: int
    sources: int
    swaps_per_edge: int
    seed: int
    delta_t: tuple[float, ...]


def initial_regular_edges(nodes: int, degree: int) -> set[tuple[int, int]]:
    """Return a simple regular circulant graph as an edge set."""
    if nodes <= degree:
        raise ValueError("nodes must be larger than degree")
    if (nodes * degree) % 2:
        raise ValueError("nodes * degree must be even")

    edges: set[tuple[int, int]] = set()
    half = degree // 2
    for u in range(nodes):
        for offset in range(1, half + 1):
            v = (u + offset) % nodes
            edges.add((min(u, v), max(u, v)))

    if degree % 2:
        if nodes % 2:
            raise ValueError("odd degree requires an even number of nodes")
        for u in range(nodes // 2):
            edges.add((u, u + nodes // 2))

    expected = nodes * degree // 2
    if len(edges) != expected:
        raise RuntimeError("failed to construct the initial regular graph")
    return edges


def random_regular_edges(
    nodes: int,
    degree: int,
    rng: np.random.Generator,
    swaps_per_edge: int,
) -> list[tuple[int, int]]:
    """Randomize a simple regular graph with degree-preserving edge swaps."""
    edge_set = initial_regular_edges(nodes, degree)
    edge_list = list(edge_set)
    target_swaps = swaps_per_edge * len(edge_list)
    accepted = 0
    attempts = 0
    max_attempts = max(100, 20 * target_swaps)

    while accepted < target_swaps and attempts < max_attempts:
        attempts += 1
        i, j = rng.choice(len(edge_list), size=2, replace=False)
        a, b = edge_list[int(i)]
        c, d = edge_list[int(j)]
        if rng.random() < 0.5:
            c, d = d, c
        if len({a, b, c, d}) < 4:
            continue

        new_1 = (min(a, d), max(a, d))
        new_2 = (min(c, b), max(c, b))
        old_1 = edge_list[int(i)]
        old_2 = edge_list[int(j)]
        if new_1 == new_2 or new_1 in edge_set or new_2 in edge_set:
            continue

        edge_set.remove(old_1)
        edge_set.remove(old_2)
        edge_set.add(new_1)
        edge_set.add(new_2)
        edge_list[int(i)] = new_1
        edge_list[int(j)] = new_2
        accepted += 1

    if accepted < target_swaps:
        raise RuntimeError(
            f"edge randomization stopped after {accepted}/{target_swaps} swaps"
        )

    degrees = np.zeros(nodes, dtype=np.int32)
    for u, v in edge_list:
        degrees[u] += 1
        degrees[v] += 1
    if not np.all(degrees == degree):
        raise RuntimeError("degree-preserving randomization failed")
    return edge_list


def poisson_events(
    edges: Iterable[tuple[int, int]],
    duration: float,
    mean_dt: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate independent stationary Poisson events for every static edge."""
    event_rows: list[tuple[float, int, int]] = []
    rate_window = duration / mean_dt
    for u, v in edges:
        count = int(rng.poisson(rate_window))
        if count:
            times = rng.uniform(0.0, duration, size=count)
            event_rows.extend((float(t), u, v) for t in times)

    event_rows.sort(key=lambda row: row[0])
    times = np.fromiter((row[0] for row in event_rows), dtype=np.float64)
    left = np.fromiter((row[1] for row in event_rows), dtype=np.int32)
    right = np.fromiter((row[2] for row in event_rows), dtype=np.int32)
    return times, left, right


def reduced_event_graph(
    times: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    nodes: int,
) -> list[list[tuple[int, float]]]:
    """Connect each event to the first later event at either endpoint."""
    adjacency: list[list[tuple[int, float]]] = [[] for _ in range(len(times))]
    last_event = np.full(nodes, -1, dtype=np.int64)

    for current, (u, v) in enumerate(zip(left, right, strict=True)):
        previous_events = {int(last_event[int(u)]), int(last_event[int(v)])}
        previous_events.discard(-1)
        for previous in previous_events:
            adjacency[previous].append(
                (current, float(times[current] - times[previous]))
            )
        last_event[int(u)] = current
        last_event[int(v)] = current

    max_out_degree = max((len(neighbours) for neighbours in adjacency), default=0)
    if max_out_degree > 2:
        raise RuntimeError("reduced event graph has out-degree larger than two")
    return adjacency


def sample_sources(
    times: np.ndarray,
    duration: float,
    count: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample sources away from both temporal boundaries."""
    candidates = np.flatnonzero(
        (times >= 0.20 * duration) & (times <= 0.40 * duration)
    )
    if not len(candidates):
        raise RuntimeError("no source events in the sampling window")
    size = min(count, len(candidates))
    return rng.choice(candidates, size=size, replace=False)


def mean_cluster_mass(
    adjacency: list[list[tuple[int, float]]],
    sources: np.ndarray,
    delta_t: float,
) -> float:
    """Return the exact mean reachable out-component size for sampled sources."""
    marks = np.zeros(len(adjacency), dtype=np.int32)
    stamp = 0
    total = 0

    for source_value in sources:
        stamp += 1
        source = int(source_value)
        marks[source] = stamp
        stack = [source]
        mass = 0
        while stack:
            event = stack.pop()
            mass += 1
            for following, waiting_time in adjacency[event]:
                if waiting_time <= delta_t and marks[following] != stamp:
                    marks[following] = stamp
                    stack.append(following)
        total += mass
    return total / len(sources)


def analytical_tau(delta_t: np.ndarray, degree: int, mean_dt: float) -> np.ndarray:
    """Equation (A5) of Badie-Modiri et al., PRE 105, 054313 (2022)."""
    k = float(degree)
    lam = float(mean_dt)
    numerator = (
        (4.0 * k**2 - 6.0 * k + 2.0) * np.exp(-k * lam * delta_t)
        + k * np.exp((1.0 - 2.0 * k) * lam * delta_t)
        - 2.0 * (k - 1.0) ** 2
    )
    return numerator / ((1.0 - 2.0 * k) * k)


def run(config: Config) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    rng = np.random.default_rng(config.seed)
    masses = np.zeros((config.realisations, len(config.delta_t)), dtype=np.float64)
    event_counts: list[int] = []
    swap_seed_values = rng.integers(0, np.iinfo(np.int64).max, size=config.realisations)

    for realization in range(config.realisations):
        local_rng = np.random.default_rng(int(swap_seed_values[realization]))
        edges = random_regular_edges(
            config.nodes, config.degree, local_rng, config.swaps_per_edge
        )
        times, left, right = poisson_events(
            edges, config.duration, config.mean_dt, local_rng
        )
        adjacency = reduced_event_graph(times, left, right, config.nodes)
        sources = sample_sources(times, config.duration, config.sources, local_rng)
        event_counts.append(len(times))

        for index, delta_t in enumerate(config.delta_t):
            masses[realization, index] = mean_cluster_mass(
                adjacency, sources, delta_t
            )

        print(
            f"realisation {realization + 1}/{config.realisations}: "
            f"{len(times)} events"
        )

    means = masses.mean(axis=0)
    if config.realisations > 1:
        errors = masses.std(axis=0, ddof=1) / math.sqrt(config.realisations)
    else:
        errors = np.zeros_like(means)

    if not np.all(np.diff(means) >= -1e-12):
        raise RuntimeError("mean cluster mass must be monotone in delta_t")
    metadata = {
        "mean_event_count": float(np.mean(event_counts)),
        "min_event_count": int(min(event_counts)),
        "max_event_count": int(max(event_counts)),
        "mass_at_smallest_delta_t": float(means[0]),
        "mass_at_largest_delta_t": float(means[-1]),
    }
    return np.asarray(config.delta_t), means, errors, metadata


def save_results(
    output_dir: Path,
    config: Config,
    delta_t: np.ndarray,
    means: np.ndarray,
    errors: np.ndarray,
    metadata: dict[str, float],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    tau = analytical_tau(delta_t, config.degree, config.mean_dt)
    theory = np.where(tau < 0.0, -1.0 / tau, np.nan)
    distance = EMPIRICAL_THRESHOLD_K9 - delta_t
    fit_mask = (delta_t >= 0.050) & (delta_t <= 0.081) & (distance > 0.0)
    fit_slope, fit_intercept = np.polyfit(
        np.log(distance[fit_mask]), np.log(means[fit_mask]), 1
    )
    effective_gamma = float(-fit_slope)

    with (output_dir / "cluster_mass.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["delta_t", "mean_cluster_mass", "standard_error", "tau", "mf_mass"]
        )
        writer.writerows(zip(delta_t, means, errors, tau, theory, strict=True))

    full_metadata = {
        "configuration": {**asdict(config), "delta_t": list(config.delta_t)},
        "validation": metadata,
        "paper_thresholds_for_k9_lambda1": {
            "analytical": ANALYTICAL_THRESHOLD_K9,
            "empirical": EMPIRICAL_THRESHOLD_K9,
        },
        "diagnostic_power_law_fit": {
            "delta_t_interval": [0.050, 0.081],
            "effective_gamma": effective_gamma,
            "warning": "diagnostic only; not a finite-size-scaling estimate",
        },
    }
    (output_dir / "run_metadata.json").write_text(
        json.dumps(full_metadata, indent=2) + "\n", encoding="utf-8"
    )

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0))

    ax = axes[0]
    ax.errorbar(
        delta_t,
        means,
        yerr=errors,
        fmt="o",
        color="#9B0014",
        markerfacecolor="white",
        capsize=2.5,
        label="simulation (mean ± SE)",
    )
    valid = np.isfinite(theory)
    ax.plot(delta_t[valid], theory[valid], color="#2764C8", label=r"mean field $-1/\tau$")
    ax.axvline(ANALYTICAL_THRESHOLD_K9, color="#2764C8", ls="--", lw=1.2,
               label=r"analytical $\delta t_c$")
    ax.axvline(EMPIRICAL_THRESHOLD_K9, color="#3A8F73", ls=":", lw=1.5,
               label=r"paper empirical $\delta t_c$")
    ax.set_xlabel(r"maximum waiting time $\delta t$")
    ax.set_ylabel(r"mean cluster mass $\langle M\rangle$")
    ax.set_yscale("log")
    ax.set_title("Subcritical growth")
    ax.legend(fontsize=8)

    ax = axes[1]
    positive = distance > 0.0
    ax.errorbar(
        distance[positive],
        means[positive],
        yerr=errors[positive],
        fmt="o",
        color="#9B0014",
        markerfacecolor="white",
        capsize=2.5,
        label="simulation",
    )
    reference_index = int(np.flatnonzero(positive)[-2]) if positive.sum() > 1 else 0
    amplitude = means[reference_index] * distance[reference_index]
    reference_x = np.geomspace(distance[positive].min(), distance[positive].max(), 200)
    ax.plot(reference_x, amplitude / reference_x, color="#3A8F73", ls="--",
            label=r"reference slope $-1$")
    fit_x = np.geomspace(distance[fit_mask].min(), distance[fit_mask].max(), 200)
    fit_y = np.exp(fit_intercept) * fit_x**fit_slope
    ax.plot(
        fit_x,
        fit_y,
        color="#2764C8",
        lw=1.6,
        label=rf"diagnostic fit $\gamma_{{\rm eff}}={effective_gamma:.2f}$",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel(r"distance from transition $\delta t_c-\delta t$")
    ax.set_ylabel(r"mean cluster mass $\langle M\rangle$")
    ax.set_title(r"Scaling toward the transition")
    ax.legend(fontsize=8)

    fig.suptitle(
        rf"Random {config.degree}-regular temporal network; Poisson activation "
        rf"($\lambda={config.mean_dt:g}$)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(output_dir / "cluster_mass.pdf", bbox_inches="tight")
    fig.savefig(output_dir / "cluster_mass.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=int, default=256)
    parser.add_argument("--degree", type=int, default=9)
    parser.add_argument("--duration", type=float, default=64.0)
    parser.add_argument("--mean-dt", type=float, default=1.0)
    parser.add_argument("--realisations", type=int, default=8)
    parser.add_argument("--sources", type=int, default=256)
    parser.add_argument("--swaps-per-edge", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument(
        "--delta-t",
        type=float,
        nargs="+",
        default=[0.020, 0.035, 0.050, 0.060, 0.068, 0.074, 0.078, 0.081, 0.083, 0.084],
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("computational/results")
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = Config(
        nodes=args.nodes,
        degree=args.degree,
        duration=args.duration,
        mean_dt=args.mean_dt,
        realisations=args.realisations,
        sources=args.sources,
        swaps_per_edge=args.swaps_per_edge,
        seed=args.seed,
        delta_t=tuple(sorted(set(args.delta_t))),
    )
    if config.degree != 9 or not math.isclose(config.mean_dt, 1.0):
        print(
            "note: threshold guide lines are literature values only for k=9, lambda=1"
        )
    delta_t, means, errors, metadata = run(config)
    save_results(args.output_dir, config, delta_t, means, errors, metadata)
    print(f"results written to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
