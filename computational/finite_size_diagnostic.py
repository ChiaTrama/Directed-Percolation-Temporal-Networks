#!/usr/bin/env python3
"""Finite-size diagnostic for the temporal-network cluster mass.

This script runs the same random 9-regular Poisson model at several spatial
sizes and compares the mean reachable event-cluster mass.  The observation
window is kept fixed, so the output is deliberately labelled a finite-size
diagnostic rather than a full finite-size-scaling analysis.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from dataclasses import asdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-complex-trovato")

import matplotlib.pyplot as plt
import numpy as np

from simulate_cluster_mass import (
    ANALYTICAL_THRESHOLD_K9,
    EMPIRICAL_THRESHOLD_K9,
    Config,
    analytical_tau,
    run,
    save_results,
)


DEFAULT_SIZES = (256, 512, 1024, 2048, 4096)
DEFAULT_DELTA_T = (0.020, 0.035, 0.050, 0.060, 0.068,
                   0.074, 0.078, 0.081, 0.083, 0.084)


def configuration_for(size: int, args: argparse.Namespace) -> Config:
    return Config(
        nodes=size,
        degree=9,
        duration=64.0,
        mean_dt=1.0,
        realisations=args.realisations,
        sources=args.sources,
        swaps_per_edge=12,
        seed=args.seed + size,
        delta_t=DEFAULT_DELTA_T,
    )


def load_dataset(result_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with (result_dir / "cluster_mass.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    delta_t = np.asarray([float(row["delta_t"]) for row in rows])
    means = np.asarray([float(row["mean_cluster_mass"]) for row in rows])
    errors = np.asarray([float(row["standard_error"]) for row in rows])
    return delta_t, means, errors


def result_matches(result_dir: Path, config: Config) -> bool:
    metadata_path = result_dir / "run_metadata.json"
    data_path = result_dir / "cluster_mass.csv"
    if not metadata_path.exists() or not data_path.exists():
        return False
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    recorded = metadata.get("configuration", {})
    expected = {**asdict(config), "delta_t": list(config.delta_t)}
    return recorded == expected


def generate_missing_sizes(
    output_dir: Path,
    sizes: tuple[int, ...],
    largest_results: Path,
    args: argparse.Namespace,
) -> dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    datasets: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    for size in sizes:
        if size == 4096:
            metadata = json.loads(
                (largest_results / "run_metadata.json").read_text(encoding="utf-8")
            )
            if int(metadata["configuration"]["nodes"]) != size:
                raise RuntimeError("the reusable largest-system result is not N=4096")
            print(f"N={size}: reusing {largest_results}", flush=True)
            datasets[size] = load_dataset(largest_results)
            continue

        config = configuration_for(size, args)
        result_dir = output_dir / f"N_{size}"
        if result_matches(result_dir, config):
            print(f"N={size}: reusing completed result", flush=True)
        else:
            print(
                f"N={size}: running {config.realisations} realizations, "
                f"{config.sources} sources",
                flush=True,
            )
            delta_t, means, errors, validation = run(config)
            save_results(result_dir, config, delta_t, means, errors, validation)
        datasets[size] = load_dataset(result_dir)

    return datasets


def save_combined_data(
    output_dir: Path,
    datasets: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> None:
    with (output_dir / "finite_size_diagnostic.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["nodes", "delta_t", "mean_cluster_mass", "standard_error"])
        for size, (delta_t, means, errors) in datasets.items():
            for values in zip(delta_t, means, errors, strict=True):
                writer.writerow([size, *values])


def save_figure(
    output_dir: Path,
    datasets: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> None:
    sizes = tuple(sorted(datasets))
    reference_delta = datasets[sizes[-1]][0]
    tau = analytical_tau(reference_delta, degree=9, mean_dt=1.0)
    theory = np.where(tau < 0.0, -1.0 / tau, np.nan)
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, len(sizes)))

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0))

    ax = axes[0]
    valid = np.isfinite(theory)
    ax.plot(
        reference_delta[valid], theory[valid], color="#2764C8", lw=2.0,
        label=r"mean field $-1/\tau$",
    )
    for color, size in zip(colors, sizes, strict=True):
        delta_t, means, errors = datasets[size]
        ax.errorbar(
            delta_t, means, yerr=errors, marker="o", ms=3.8, lw=1.2,
            capsize=1.8, color=color, label=rf"$N={size}$",
        )
    ax.axvline(
        ANALYTICAL_THRESHOLD_K9, color="#2764C8", ls="--", lw=1.1,
        label=r"analytical $\delta t_c$",
    )
    ax.axvline(
        EMPIRICAL_THRESHOLD_K9, color="#3A8F73", ls=":", lw=1.4,
        label=r"paper empirical $\delta t_c$",
    )
    ax.set_yscale("log")
    ax.set_xlabel(r"maximum waiting time $\delta t$")
    ax.set_ylabel(r"mean cluster mass $\langle M\rangle$")
    ax.set_title("Mass growth for different system sizes")
    ax.legend(fontsize=7.2, ncol=2)

    ax = axes[1]
    largest_means = datasets[sizes[-1]][1]
    largest_errors = datasets[sizes[-1]][2]
    selected_delta = (0.074, 0.081, 0.084)
    selected_colors = ("#2764C8", "#E58A32", "#9B0014")
    for target, color in zip(selected_delta, selected_colors, strict=True):
        index = int(np.argmin(np.abs(reference_delta - target)))
        ratios = []
        ratio_errors = []
        for size in sizes:
            _, means, errors = datasets[size]
            ratio = means[index] / largest_means[index]
            ratios.append(ratio)
            if size == sizes[-1]:
                ratio_errors.append(0.0)
            else:
                relative_variance = (
                    (errors[index] / means[index]) ** 2
                    + (largest_errors[index] / largest_means[index]) ** 2
                )
                ratio_errors.append(ratio * np.sqrt(relative_variance))
        ax.errorbar(
            sizes, ratios, yerr=ratio_errors, marker="o", capsize=2.5,
            color=color, lw=1.5, label=rf"$\delta t={target:.3f}$",
        )
    ax.axhline(1.0, color="#555D68", ls="--", lw=1.1,
               label=r"$N=4096$ reference")
    ax.set_xscale("log", base=2)
    ax.set_xticks(sizes, [str(size) for size in sizes])
    ax.set_xlabel(r"number of nodes $N$")
    ax.set_ylabel(r"mass ratio $\langle M\rangle_N/\langle M\rangle_{4096}$")
    ax.set_title("Finite-size suppression near the transition")
    ax.legend(fontsize=8)

    fig.suptitle(
        r"Finite-size diagnostic at fixed $T=64$: random 9-regular, Poisson "
        r"activation ($\lambda=1$)",
        fontsize=11.5,
    )
    fig.tight_layout()
    fig.savefig(output_dir / "finite_size_diagnostic.pdf", bbox_inches="tight")
    fig.savefig(output_dir / "finite_size_diagnostic.png", dpi=220,
                bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--realisations", type=int, default=64)
    parser.add_argument("--sources", type=int, default=1024)
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("computational/results_finite_size"),
    )
    parser.add_argument(
        "--largest-results", type=Path,
        default=Path("computational/results_high_statistics"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    datasets = generate_missing_sizes(
        args.output_dir, DEFAULT_SIZES, args.largest_results, args
    )
    save_combined_data(args.output_dir, datasets)
    save_figure(args.output_dir, datasets)
    metadata = {
        "sizes": list(DEFAULT_SIZES),
        "fixed_duration": 64.0,
        "realisations_for_N_256_to_2048": args.realisations,
        "sources_per_realisation": args.sources,
        "largest_result_reused_from": str(args.largest_results),
        "interpretation": "finite-size diagnostic, not full finite-size scaling",
    }
    (args.output_dir / "finite_size_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"combined results written to {args.output_dir.resolve()}", flush=True)


if __name__ == "__main__":
    main()
