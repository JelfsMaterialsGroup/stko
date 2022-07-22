import re
from os.path import join

from setuptools import find_packages, setup


def get_version():
    with open(join('stko', '__init__.py'), 'r') as f:
        content = f.read()
    p = re.compile(r'^__version__ = [\'"]([^\'\"]*)[\'"]', re.M)
    return p.search(content).group(1)


setup(
    name="stko",
    version="0.0.40",
    author="Steven Bennett, Andrew Tarzia",
    description=(
        "Contains molecular optimizers and property calculators for "
        "use with stk."
    ),
    url="https://github.com/JelfsMaterialsGroup/stko",
    packages=find_packages(),
    install_requires=(
        'scipy',
        'numpy',
        'networkx',
        'stk',
        'rdkit-pypi',
    ),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
)
