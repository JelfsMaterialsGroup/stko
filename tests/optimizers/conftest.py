import pytest
import numpy as np
import stko
import stk


def a_molecule():
    return stk.BuildingBlock(smiles='CC')


class PassingOptimizer(stko.Optimizer):

    def optimize(self, mol):
        return a_molecule().with_centroid(np.array(([1, 3, 3])))


class FailingOptimizer(stko.Optimizer):

    def optimize(self, mol):
        raise Exception()


@pytest.fixture
def passing_optimizer():
    return PassingOptimizer()


@pytest.fixture
def failing_optimizer():
    return FailingOptimizer()


class CaseData:
    """
    A test case.

    Attributes:
        molecule:
            The molecule to be tested.

        unoptimised_energy:
            The energy of the molecule from stk generation.

        optimised_energy:
            The energy of the molecule after optimisation.

    """

    position_matrix: np.ndarray

    def __init__(
        self,
        molecule,
        unoptimised_energy,
        optimised_energy,
    ):

        self.molecule = molecule
        self.unoptimised_energy = unoptimised_energy
        self.optimised_energy = optimised_energy


@pytest.fixture(
    scope='session',
    params=[
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            unoptimised_energy=98.68403256936926,
            optimised_energy=43.16966811699349,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)'
                'cc3c2c1'
            ),
            unoptimised_energy=2032.8111567743895,
            optimised_energy=805.7039129507456,
        ),
        CaseData(
            molecule=stk.BuildingBlock('CCCCCC'),
            unoptimised_energy=141.44622279628743,
            optimised_energy=20.692758531646255,
        ),
        CaseData(
            molecule=stk.BuildingBlock('c1ccccc1'),
            unoptimised_energy=56.73128282588534,
            optimised_energy=44.27000519453281,
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCNCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='AB',
                    num_repeating_units=2,
                    optimizer=stk.MCHammer(),
                ),
            ),
            unoptimised_energy=15002.946293854524,
            optimised_energy=110.38365695234553,
        ),
    ],
)
def case_molecule(request):
    """
    A :class:`.Molecule` instance.

    """

    return request.param
