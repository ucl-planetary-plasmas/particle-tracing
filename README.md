# Particle tracing

This is a software package to study the motion of charged particles in
the magnetosphere of Jupiter and Saturn using the AGA magnetodisc model

This repository contains both:

- a MATLAB implementation
- a Python implementation


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

## AGA Magnetodisc models installation

Run the commands

```bash
cd data
./downloadMDiscFiles.sh
```

## MATLAB version

## Python version
