
# Particle tracing framework (MATLAB + Python)

This repository contains two independent implementations (MATLAB and Python)
for simulating and analysing motion of charged particles in the gas giants
(Jupiter and Saturn) magetosphere using the AGA magnetodisc model, along
with shared datasets and utiliies

## Overview

The project provides tools to:

- Visualise magnetic field models for planetary magnetodiscs
- Trace charged particle trajectories
- Compute bounce and drift periods
- Compare results between MATLAB and Python implementations
- Run example simulations for Earth, Jupiter, and Saturn

---

## Repository Structure

```text
particle-tracing/
├── data/
│   └── downloadMDiscFiles.sh
│
├── matlab/
│   ├── cmptraj_earth.m
│   ├── cmptraj_jupiter.m
│   ├── cmptraj_saturn.m
│   ├── dipbtracer.m
│   ├── dipoleMagneticField3D.m
│   ├── diptracer.m
│   ├── example_trajectory.m
│   ├── fitTb.m
│   ├── fitTd.m
│   ├── getBmoddipole.m
│   ├── getbounceperiod.m
│   ├── getBvectordipole.m
│   ├── getdriftperiod.m
│   ├── leasqr.m
│   ├── mdbtracer.m
│   ├── MDiscField.m
│   ├── mdiscMagneticField3D.m
│   ├── mdtracer.m
│   ├── runs.m
│   ├── testmdisc.m
│   └── trajectory_main.m
│
├── python/
│   ├── dist/
│   ├── examples/
│   ├── LICENSE.txt
│   ├── notebooks/
│   ├── pymagdisc/
│   ├── pyproject.toml
│   ├── README.md
│   ├── requirements.txt
│   ├── setup.cfg
│   ├── setup.py
│   └── tests/
│
└── README.md
```

---

## Data

The `data/` directory contains scripts for downloading required magnetodisc
datasets (AGA Magnetodisc models)

```bash
bash data/downloadMDiscFiles.sh
```

---

## MATLAB Usage


### Run example simulations

```bash
cd matlab
matlab -batch "example_trajectory"
```

### Run full trajectory simulation

```bash
matlab -batch "trajectory_main"
```

---

## Python Usage


### Installation

```bash
cd python
pip install .
```

### Development installation

```bash
pip install -e .
```

### Run tests

```bash
pytest
```

### Run examples

```bash
python examples/example.py
```

---

## Requirements

### MATLAB
- MATLAB R20XX+
- Required toolboxes (if applicable)

### Python
- Python 3.7+
- pip

---

## Features

- Multi-planet magnetodisc modelling (Earth, Jupiter, Saturn)
- Particle trajectory tracing
- Analytical and numerical field models
- Cross-validation between MATLAB and Python implementations

---

## Citation

If you use this software, please cite:

> Add paper / DOI / reference here

---

## License

Add license information here.

---

## Authors

- Patrick Guio
- Contributors

---

## Acknowledgements


