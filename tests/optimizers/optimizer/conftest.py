import numpy as np
import pytest
import stk
import stko


def a_molecule():
    return stk.BuildingBlock(smiles="CC")


@pytest.fixture
def unoptimized_mol():
    return a_molecule()


class PassingOptimizer(stko.Optimizer):
    def optimize(self, mol: stk.Molecule) -> stk.Molecule:
        return a_molecule().with_centroid(np.array(([1, 3, 3])))


class FailingOptimizer(stko.Optimizer):
    def optimize(self, mol: stk.Molecule) -> stk.Molecule:
        raise Exception()


@pytest.fixture
def passing_optimizer():
    return PassingOptimizer()


@pytest.fixture
def failing_optimizer():
    return FailingOptimizer()
