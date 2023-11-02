# List all recipes.
default:
  @just --list

# Build docs.
docs:
  rm -rf docs/source/_autosummary
  make -C docs html
  echo Docs are in $PWD/docs/build/html/index.html

# Install development environment.
dev:
  pip install -e '.[dev]'

# Run code checks.
check:
  #!/usr/bin/env bash

  error=0
  trap error=1 ERR

  echo
  (set -x; ruff . )

  echo
  ( set -x; black --check . )

  echo
  ( set -x; mypy src )

  echo
  ( set -x; pytest --cov=stko --cov-report term-missing )

  test $error = 0

# Auto-fix code issues.
fix:
  black .
  ruff --fix .
