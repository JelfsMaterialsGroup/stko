import numpy as np
import pytest
import stk

import stko


def a_molecule() -> stk.BuildingBlock:
    return stk.BuildingBlock(smiles="CC")


@pytest.fixture
def unoptimized_mol() -> stk.BuildingBlock:
    return a_molecule()


class PassingOptimizer(stko.Optimizer):
    def optimize(self, mol: stko.MoleculeT) -> stko.MoleculeT:
        return mol.with_centroid(np.array([1, 3, 3]))


class FailingOptimizer(stko.Optimizer):
    def optimize(self, mol: stko.MoleculeT) -> stko.MoleculeT:  # noqa: ARG002
        raise RuntimeError


@pytest.fixture
def passing_optimizer() -> PassingOptimizer:
    return PassingOptimizer()


@pytest.fixture
def failing_optimizer() -> FailingOptimizer:
    return FailingOptimizer()
