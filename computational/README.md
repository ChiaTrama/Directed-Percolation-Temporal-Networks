# Computational illustration

## Goal

Reproduce, at a deliberately smaller scale, the qualitative result shown in
Fig. 7(b) of Badie-Modiri et al., *Phys. Rev. E* **105**, 054313 (2022): the
growth of the mean reachable event-cluster mass as the limited waiting time
approaches the transition from the subcritical side.

This is an **illustration**, not a publication-scale reproduction. The paper
uses systems up to `N = 2^17` and hundreds or thousands of realizations. The
defaults below are intended to run on a laptop and to produce a result that can
be understood and explained during the exam.

## Model

1. Construct an approximately uniform random `k`-regular static network. The
   implementation starts from a regular circulant graph and randomizes it with
   degree-preserving double-edge swaps.
2. Activate every static edge with an independent Poisson process of mean
   inter-event time `lambda = 1`.
3. Sort all temporal events by time.
4. Construct the reduced event graph: from each event, retain only the first
   later event involving each of its two endpoints.
5. For every tested waiting time `delta_t`, retain reduced-event-graph edges
   whose waiting time is at most `delta_t`.
6. From randomly sampled source events, measure the number of reachable events
   `M` and average it over sources and network realizations.

The analytical control parameter for a random `k`-regular base network with
Poisson activation is Eq. (A5) of the PRE paper:

```text
tau(delta_t) = [
    (4 k^2 - 6 k + 2) exp(-k lambda delta_t)
    + k exp((1 - 2 k) lambda delta_t)
    - 2 (k - 1)^2
] / [(1 - 2 k) k]
```

In the mean-field subcritical regime, `M = -1/tau`, giving `gamma = 1`.

For `k = 9` and `lambda = 1`, the paper reports:

- analytical threshold: `delta_t_c = 0.08559`;
- empirical threshold: `delta_t_c = 0.08808`.

## Run

From the project root:

```bash
python3 computational/simulate_cluster_mass.py
```

Faster smoke test (writes to a separate output directory):

```bash
python3 computational/simulate_cluster_mass.py \
  --nodes 128 --duration 32 --realisations 2 --sources 64 \
  --output-dir computational/results_smoke
```

Outputs are written to `computational/results/`:

- `cluster_mass.csv`: numerical values and uncertainty;
- `cluster_mass.png`: presentation-ready raster figure;
- `cluster_mass.pdf`: vector version for Beamer;
- `run_metadata.json`: parameters and basic validation information.

## High-statistics run used in the presentation

The figure included in the revised presentation was generated with:

```bash
python3 computational/simulate_cluster_mass.py \
  --nodes 4096 \
  --realisations 128 \
  --sources 1024 \
  --output-dir computational/results_high_statistics
```

This run processes approximately 1.18 million temporal events per realization.
On the development machine it completed in about 22 minutes and 38 seconds.
The original exploratory run is omitted from this release; `computational/results/`
is created when you run the default simulation.

## Finite-size diagnostic

To compare the same observable across several spatial system sizes, run:

```bash
python3 -u computational/finite_size_diagnostic.py
```

The script simulates `N = 256, 512, 1024, 2048` with 64 realizations and 1024
sampled sources, and reuses the high-statistics `N = 4096` result. Combined
outputs are written to `computational/results_finite_size/`. The observation
duration is fixed at `T = 64`, so this comparison demonstrates finite-size
suppression but is not a full finite-size-scaling or data-collapse analysis.

## Interpretation and limitations

- The relevant result is the rapid growth of `M` when approaching the critical
  region from below.
- Agreement with `M = -1/tau` is expected to be qualitative because the
  analytical theory treats the event graph as a random directed graph with the
  same degree distribution.
- Deviations near the transition are expected from finite network size, finite
  observation time, temporal-boundary effects, and the small number of
  realizations.
- This experiment alone is not enough to claim a precise measurement of
  `gamma`. A robust exponent estimate would require larger sizes, many more
  realizations and a finite-size-scaling analysis.
- The plot reports a diagnostic effective exponent fitted over
  `0.050 <= delta_t <= 0.081`. Its purpose is to check consistency with
  `gamma = 1`, not to provide an independent precision measurement.
