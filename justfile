# List all recipes.
default:
  @just --list

# Build docs.
docs:
  rm -rf ./docs/build docs/source/_autosummary
  make -C docs html
  echo Docs are in $PWD/docs/build/html/index.html

# Install development environment.
dev:
  pip install -e '.[dev]'
  mamba install -y -c conda-forge pytorch==2.3.1 torchdata==0.7.1 openff-toolkit openmm openmmtools rdkit!=2024.3.5 dgl espaloma_charge
  mamba install -y xtb
  mamba install -y openbabel
  mamba install -y mdanalysis


# Run code checks.
check:
  #!/usr/bin/env bash

  error=0
  trap error=1 ERR

  echo
  (set -x; ruff check . )

  echo
  ( set -x; ruff format --check . )

  echo
  ( set -x; mypy src examples )

  echo
  ( set -x; pytest --cov=stko --cov-report term-missing )

  echo
  ( set -x; make -C docs doctest )

  test $error = 0

# Auto-fix code issues.
fix:
  ruff format .
  ruff check --fix .
