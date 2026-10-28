# Directed Percolation in Temporal Networks

**Physics of Complex Systems · MSc in Physics of Data · University of Padua · September 2026**<br>
Chiara Tramarin

Course project based on the study and presentation of two papers on directed
percolation in temporal networks. The presentation reviews temporal paths,
event graphs, absorbing-state transitions and the mean-field description of
critical behaviour.

A small Python simulation accompanies the presentation and qualitatively
illustrates the growth of the mean reachable event-cluster mass near the
transition. It is intended as a course exercise rather than a complete
reproduction of the original numerical study.

[View the final presentation](presentation/directed_percolation.pdf) ·
[Explore the computational illustration](computational/README.md)

## Quick start

Use **Python 3.10 or later**. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python computational/simulate_cluster_mass.py \
  --nodes 128 --duration 32 --realisations 2 --sources 64 \
  --output-dir computational/results_smoke
```

On Windows, activate the environment with `.venv\Scripts\activate` instead.
This is a small workflow check. Full commands, saved datasets and modelling
limitations are documented in the
[computational README](computational/README.md).

## Repository contents

- `presentation/directed_percolation.pdf`: final presentation.
- `computational/`: simulation scripts, documentation and saved numerical results.

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
