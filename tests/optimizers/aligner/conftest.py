import numpy as np
import pytest
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
    scope="session",
    params=[
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            initial_molecule=(
                stk.BuildingBlock("NCCN").with_rotation_about_axis(
                    1.34,
                    np.array((0, 0, 1)),
                    np.array((0, 0, 0)),
                )
            ),
            rmsd=0.22366328852274148,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            initial_molecule=(
                stk.BuildingBlock("CCCCCC").with_rotation_about_axis(
                    1.34,
                    np.array((0, 0, 1)),
                    np.array((0, 0, 0)),
                )
            ),
            rmsd=0.12594321424488517,
        ),
    ],
)
def case_molecule(request):
    """
    A :class:`stk.Molecule` instance.

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
    scope="session",
    params=[
        CasePotential(
            molecule=stk.BuildingBlock("NCCN"),
            initial_molecule=(
                stk.BuildingBlock("NCCN").with_rotation_about_axis(
                    1.34,
                    np.array((0, 0, 1)),
                    np.array((0, 0, 0)),
                )
            ),
            potential=16.33559420081716,
            pairs=(("C", "C"), ("N", "N")),
        ),
        CasePotential(
            molecule=stk.BuildingBlock("NCCN"),
            initial_molecule=(
                stk.BuildingBlock("NCCN").with_rotation_about_axis(
                    1.34,
                    np.array((0, 0, 1)),
                    np.array((0, 0, 0)),
                )
            ),
            potential=4.037489677500126,
            pairs=(("C", "C"),),
        ),
        CasePotential(
            molecule=stk.BuildingBlock("NCCN"),
            initial_molecule=(
                stk.BuildingBlock("NCCN").with_rotation_about_axis(
                    1.34,
                    np.array((0, 0, 1)),
                    np.array((0, 0, 0)),
                )
            ),
            potential=12.298104523317033,
            pairs=(("N", "N"),),
        ),
        CasePotential(
            molecule=stk.BuildingBlock("NCCN"),
            initial_molecule=(
                stk.BuildingBlock("NCCN").with_rotation_about_axis(
                    1.34,
                    np.array((0, 0, 1)),
                    np.array((0, 0, 0)),
                )
            ),
            potential=0.0,
            pairs=(),
        ),
    ],
)
def case_potential(request):
    """
    A :class:`stk.Molecule` instance.

    """

    return request.param
