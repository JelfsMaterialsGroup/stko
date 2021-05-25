import pytest
import numpy as np
import stko
import stk


class PassingOptimizer(stko.Optimizer):

    def optimize(self, mol):
        atoms = [stk.Atom(0, 1), stk.Atom(1, 1)]
        return stk.Molecule(
            atoms=atoms,
            bonds=[stk.Bond(atoms[0], atoms[1], 1)],
            position_matrix=np.array(([0, 0, 0], [0, 1, 0])),
        )


class FailingOptimizer(stko.Optimizer):

    def optimize(self, mol):
        raise Exception()


@pytest.fixture
def passing_optimizer():
    return PassingOptimizer()


@pytest.fixture
def failing_optimizer():
    return FailingOptimizer()


@pytest.fixture
def unoptimized_mol():
    return stk.Molecule(
        atoms=[stk.C(0, 1)],
        bonds=[],
        position_matrix=np.array([0, 0, 0]),
    )


@pytest.fixture
def optimized_mol():
    atoms = [stk.Atom(0, 1), stk.Atom(1, 1)]
    return stk.Molecule(
        atoms=atoms,
        bonds=[stk.Bond(atoms[0], atoms[1], 1)],
        position_matrix=np.array(([0, 0, 0], [0, 1, 0])),
    )
