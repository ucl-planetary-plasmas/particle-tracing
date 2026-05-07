from pathlib import Path
from setuptools import setup, find_packages

requirements = Path("requirements.txt").read_text().splitlines()

install_requires = [
    line.strip()
    for line in requirements
    if line.strip() and not line.startswith("#")
]

setup(
    name="pymagdisc",
    version="1.0.0",
    description="Python package for magnetodisc model.",
    author="I Kit Cheng, Dimitrios Millas",
    author_email="i.cheng.19@ucl.ac.uk, dimitrios.millas@ucl.ac.uk",

    packages=find_packages(exclude=["tests", "tests.*"]),

    install_requires=install_requires,

    entry_points={
        "console_scripts": [
            "pymagdisc-plot_mdisc=pymagdisc.vis.plot_magdisc:plot_command",
        ]
    },

    python_requires=">=3.9",
)
