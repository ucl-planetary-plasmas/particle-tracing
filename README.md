
# Charged Particle Tracing in Planetary Magnetospheres

[![MATLAB](https://img.shields.io/badge/MATLAB-R20XX+-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

The repository `particle-tracing` contains two independent implementations
(MATLAB and Python) for simulating and analysing charged particle motion in
planetary magnetospheres. 

Available planets are the gas giants (Jupiter and
Saturn) with the AGA (Achilleos, Guio and Arridge, 2010; [link](https://doi.org/10.1111/j.1365-2966.2009.15865.x)) magnetic field model (magnetodisc),
and the Earth with a simple dipole magnetic field. 

`particle-tracing` includes magnetic field models, particle tracing tools,
and analysis utilities.


## 🚀 Key Features

The project provides tools to:

- Visualise magnetic field models for planetary magnetodiscs
- Trace charged particle trajectories
- Compute bounce and drift periods
- Compare results between MATLAB and Python implementations
- Run example simulations for Earth, Jupiter, and Saturn

---

## 📁 Repository Structure

```text
particle-tracing/
├── data/                 # Dataset download utilities
│   └── downloadMDiscFiles.sh
│
├── matlab/              # MATLAB implementation
│   ├── cmptraj_earth.m
│   ├── cmptraj_jupiter.m
│   ├── cmptraj_saturn.m
│   ├── dipbtracer.m
│   ├── mdbtracer.m
│   ├── runs.m
│   ├── testmdisc.m
│   ├── example_trajectory.m
│   └── ...

│
├── python/              # Python implementation
│   ├── pymagdisc/
│   ├── tests/
│   ├── examples/
│   ├── notebooks/
│   ├── pyproject.toml
│   ├── setup.py
│   └── ...
│
└── README.md
```

---

## 📦 Data Setup

Download required magnetodisc datasets: the `data/` directory contains
scripts for downloading required magnetodisc datasets (AGA Magnetodisc
models)

```bash id="data_setup"
bash data/downloadMDiscFiles.sh
```

---

## 🧪 MATLAB Usage


### Run example simulations

```bash id="matlab_example"
cd matlab
matlab -batch "example_trajectory"
```


---

## 🐍 Python Usage


### Installation

```bash id="py_install"
cd python
pip install .
```

### Development install

```bash id="py_dev"
pip install -e .
```

### Run tests

```bash id="py_tests"
pytest
```

### Run examples

```bash id="py_examples"
python examples/example_trajectory.py
```

---

## 📌 Requirements

### MATLAB
- MATLAB R20XX or newer
- Required toolboxes (if applicable)

### Python
- Python 3.7+
- pip

---

## 📊 Scientific Capabilities

- Multi-planet magnetic field modelling (Earth, Jupiter, Saturn)
- Particle trajectory tracing
- Analysis / diagnostics of particle tracing

---

## 📖 Citation

If you use this software, please cite:

> Add paper / DOI / reference here

---

## 🤝 Contributing

Contributions are welcome.

1. Fork repository
2. Create feature branch
3. Commit changes
4. Open pull request

---

## 📄 License

Add license information here (e.g., MIT, GPL, etc.).

---

## 👤 Authors

- Author: Patrick Guio <p.guio@ucl.ac.uk>
- Contributors: 
    - I Kit Cheng <i.cheng.19@ucl.ac.uk> 
    - Dimitrios Millas <dimitrios.millas@ucl.ac.uk>



