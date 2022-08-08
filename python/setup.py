from setuptools import setup, find_packages

with open('./requirements.txt', 'r') as f:
    reqs = f.read().splitlines()
    install_requires = [req.split('=')[0] for req in reqs]

setup(
    name="pymagdisc",
    author="I Kit Cheng, Dimitrios Millas",
    author_email="i.cheng.19@ucl.ac.uk, dimitrios.millas@ucl.ac.uk",
    description="Python package for magnetodisc model.",
    version="1.0",
    packages=find_packages(exclude=["tests.*", "tests"]),
    install_requires=install_requires,
    setup_requires=(
        'pytest-runner',
        ),
    tests_require=(
            'pytest-cov',
        ),

    entry_points={
       'console_scripts': [
           'pymagdisc-plot_mdisc = pymagdisc.vis.plot_magdisc:plot_command',
       ]
    },
)