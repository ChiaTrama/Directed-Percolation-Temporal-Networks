# Directed Percolation in Temporal Networks

**Physics of Complex Systems · MSc in Physics of Data · University of Padua · September 2026**<br>
Chiara Tramarin

This repository contains the materials prepared for a Physics of Complex Systems
exam project on directed percolation in temporal networks. The main outcome is a
literature-based presentation of the papers cited below, accompanied by a small
Python computational illustration of one qualitative result.

[View the final presentation](presentation/directed_percolation.pdf) ·
[Computational methods and reproduction commands](computational/README.md)

## Project scope

The presentation reviews temporal paths, event graphs, absorbing-state
transitions, mean-field theory and critical behaviour. Its goal is to explain
how a limited waiting time maps reachability in temporal networks to a directed
percolation problem.

The theory and reference thresholds come from the cited papers. My contribution
was the presentation and a small-scale Python implementation used to illustrate
the reported growth of the mean reachable event-cluster mass near the transition.
The numerical part is a qualitative course exercise rather than an independent
research result or a complete reproduction of the original study.

**Methods and tools:** Python, NumPy, Matplotlib, random graph generation,
stochastic simulation, graph traversal, uncertainty estimation and scientific
visualisation.

## Computational illustration

The accompanying code generates random regular networks with Poisson edge
activations and measures the mean number of events reachable under a maximum
waiting-time constraint. Saved data and figures include a higher-statistics run
and a comparison across several network sizes. Parameters, random seeds, results
and limitations are documented in the
[computational README](computational/README.md).

## Quick start

Use **Python 3.10 or later**. The dependencies are pinned to the versions used
for the local packaging checks (Python 3.10.12).

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python computational/simulate_cluster_mass.py \
  --nodes 128 --duration 32 --realisations 2 --sources 64 \
  --output-dir computational/results_smoke
```

On Windows, activate the environment with `.venv\Scripts\activate` instead.
The simulation writes CSV data, PNG/PDF plots and JSON metadata. This small run
checks the workflow; it does not reproduce the statistics of the saved figures.

To regenerate the combined finite-size plot from the included datasets:

```bash
python computational/finite_size_diagnostic.py
```

With the default parameters, the existing datasets are reused. See the
[computational README](computational/README.md) for the longer simulation command
used to generate the presentation result and for modelling limitations.

## Repository contents

- `presentation/directed_percolation.pdf`: final presentation (v4).
- `computational/simulate_cluster_mass.py`: network generation and cluster-mass simulation.
- `computational/finite_size_diagnostic.py`: system-size comparison.
- `computational/results_high_statistics/`: presentation figure, data and metadata.
- `computational/results_finite_size/`: combined figure and numerical inputs for each size.

## References

1. A. Badie-Modiri, A. K. Rizi, M. Karsai and M. Kivelä,
   *Directed percolation in temporal networks*, Physical Review Research **4**,
   L022047 (2022). [Paper](https://doi.org/10.1103/PhysRevResearch.4.L022047).
2. A. Badie-Modiri, A. K. Rizi, M. Karsai and M. Kivelä,
   *Directed percolation in random temporal network models with heterogeneities*,
   Physical Review E **105**, 054313 (2022).
   [Paper](https://doi.org/10.1103/PhysRevE.105.054313).

The presentation includes discussion of these works and attributed reference
material. Paper PDFs are not redistributed in this repository.
