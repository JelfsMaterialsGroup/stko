import pytest
import numpy as np
import stk


class CaseData:
    """
    A test case.

    Attributes:
        molecule:
            The molecule to be tested.

        initial_molecule:
            The initial molecule to be tested.

        rmsd:
            The rmsd between two molecules.

    """

    position_matrix: np.ndarray

    def __init__(self, molecule, initial_molecule, rmsd):

        self.molecule = molecule
        self.initial_molecule = initial_molecule
        self.rmsd = rmsd


@pytest.fixture(
    scope='session',
    params=[
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            initial_molecule=(
                stk.BuildingBlock('NCCN').with_rotation_about_axis(
                    1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
                )
            ),
            rmsd=0.2440040878642247,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)'
                'cc3c2c1'
            ),
            initial_molecule=(
                stk.BuildingBlock(
                    'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cncc'
                    'c54)cc3c2c1'
                ).with_rotation_about_axis(
                    1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
                )
            ),
            rmsd=0.06634296595541135,
        ),
        CaseData(
            molecule=stk.BuildingBlock('CCCCCC'),
            initial_molecule=(
                stk.BuildingBlock('CCCCCC').with_rotation_about_axis(
                    1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
                )
            ),
            rmsd=0.21668989505225564,
        ),
        CaseData(
            molecule=stk.BuildingBlock('c1ccccc1'),
            initial_molecule=(
                stk.BuildingBlock('c1ccccc1').with_rotation_about_axis(
                    1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
                )
            ),
            rmsd=0.10509500676442843,
        ),
    ],
)
def case_molecule(request):
    """
    A :class:`.Molecule` instance.

    """

    return request.param


class CasePotential:
    """
    A test case.

    Attributes:
        molecule:
            The molecule to be tested.

        initial_molecule:
            The initial molecule to be tested.

        potential:
            The calculated potential.

        pairs:
            The atom pair types to calculate potential between.

    """

    position_matrix: np.ndarray

    def __init__(self, molecule, initial_molecule, potential, pairs):

        self.molecule = molecule
        self.initial_molecule = initial_molecule
        self.potential = potential
        self.pairs = pairs


@pytest.fixture(
    scope='session',
    params=[
        CasePotential(
            molecule=stk.BuildingBlock('NCCN'),
            initial_molecule=(
                stk.BuildingBlock('NCCN').with_rotation_about_axis(
                    1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
                )
            ),
            potential=16.33559420081716,
            pairs=(('C', 'C'), ('N', 'N')),
        ),
        CasePotential(
            molecule=stk.BuildingBlock('NCCN'),
            initial_molecule=(
                stk.BuildingBlock('NCCN').with_rotation_about_axis(
                    1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
                )
            ),
            potential=4.037489677500126,
            pairs=(('C', 'C'), ),
        ),
        CasePotential(
            molecule=stk.BuildingBlock('NCCN'),
            initial_molecule=(
                stk.BuildingBlock('NCCN').with_rotation_about_axis(
                    1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
                )
            ),
            potential=12.298104523317033,
            pairs=(('N', 'N'), ),
        ),
        CasePotential(
            molecule=stk.BuildingBlock('NCCN'),
            initial_molecule=(
                stk.BuildingBlock('NCCN').with_rotation_about_axis(
                    1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
                )
            ),
            potential=0.0,
            pairs=(),
        ),
    ],
)
def case_potential(request):
    """
    A :class:`.Molecule` instance.

    """

    return request.param
