from dataclasses import dataclass

import numpy as np
import pytest
import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    initial_molecule: stk.Molecule
    rmsd: float


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
            rmsd=0.1996375974321242,
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
            rmsd=0.14932987908894566,
        ),
    ],
)
def case_molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param


@dataclass(frozen=True, slots=True)
class CasePotential:
    molecule: stk.Molecule
    initial_molecule: stk.Molecule
    potential: float
    pairs: tuple[tuple[str, str], ...]


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
            potential=16.937605187959882,
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
            potential=3.897053132727992,
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
            potential=13.04055205523189,
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
def case_potential(request: pytest.FixtureRequest) -> CasePotential:
    return request.param
