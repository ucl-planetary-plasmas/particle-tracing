from pathlib import Path
from setuptools import setup, find_packages

README = Path("README.md").read_text(encoding="utf-8")

setup(
    name="pymagdisc",
    version="1.0.0",
    description="Python package for magnetodisc model.",
    long_description=README,
    long_description_content_type="text/markdown",

    author="I Kit Cheng, Dimitrios Millas",
    author_email="i.cheng.19@ucl.ac.uk, dimitrios.millas@ucl.ac.uk",

    packages=find_packages(exclude=["tests", "tests.*"]),

    install_requires=[
        "numpy>=1.22",
        "scipy>=1.10",
        "matplotlib>=3.3",
        "tqdm>=4.66",
    ],

    entry_points={
        "console_scripts": [
            "pymagdisc-plot_mdisc=pymagdisc.vis.plot_magdisc:plot_command",
        ]
    },

    python_requires=">=3.9",
)
