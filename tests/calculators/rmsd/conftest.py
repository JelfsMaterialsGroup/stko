from dataclasses import dataclass

import numpy as np
import pytest
import stk

import stko


@dataclass(frozen=True, slots=True)
class CaseData:
    mol1: stk.Molecule
    mol2: stk.Molecule
    rmsd: float
    kabsch_rmsd: float


_optimizer = stko.UFF()


_cc_molecule = stk.BuildingBlock("[C][C]")


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_centroid(np.array((4, 0, 0))),
            rmsd=0.0,
            kabsch_rmsd=0.0,
        ),
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_position_matrix(
                np.array(
                    [[0.7520009, 0.0, 0.0], [-0.7520009, 0.0, 0.0]],
                )
            ),
            rmsd=0.0,
            kabsch_rmsd=0.0,
        ),
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_position_matrix(
                np.array(
                    [[1.7520009, 0.0, 0.0], [-1.7520009, 0.0, 0.0]],
                )
            ),
            rmsd=1.0,
            kabsch_rmsd=1.0,
        ),
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=_optimizer.optimize(stk.BuildingBlock("NCCN")),
            rmsd=0.191569942286972,
            kabsch_rmsd=0.17921849894135589,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=stk.BuildingBlock("CCCCCC"),
            rmsd=0.0,
            kabsch_rmsd=0.0,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=_optimizer.optimize(stk.BuildingBlock("CCCCCC")),
            rmsd=0.5178064489458448,
            kabsch_rmsd=0.42179363084243615,
        ),
        CaseData(
            mol1=stk.BuildingBlock("c1ccccc1"),
            mol2=_optimizer.optimize(stk.BuildingBlock("c1ccccc1")),
            rmsd=0.023048649062610855,
            kabsch_rmsd=0.02304785922542081,
        ),
    ],
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    """A pair of :class:`stk.Molecule` instances and an RMSD."""
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=_optimizer.optimize(stk.BuildingBlock("NCCN")),
            rmsd=0.11778216179942502,
            kabsch_rmsd=0.20811702035676308,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=_optimizer.optimize(stk.BuildingBlock("CCCCCC")),
            rmsd=0.3219861633412275,
            kabsch_rmsd=0.22563756374632568,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=stk.BuildingBlock("CCCCCC"),
            rmsd=0.0,
            kabsch_rmsd=0.0,
        ),
        CaseData(
            mol1=stk.BuildingBlock("c1ccccc1"),
            mol2=_optimizer.optimize(stk.BuildingBlock("c1ccccc1")),
            rmsd=0.021116969610025644,
            kabsch_rmsd=0.029156836455717483,
        ),
    ],
)
def ignore_h_case_data(request: pytest.FixtureRequest) -> CaseData:
    """A pair of :class:`stk.Molecule` instances and an RMSD."""
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=stk.BuildingBlock("CCCCCC"),
            rmsd=0.0,
            kabsch_rmsd=0.0,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=stk.BuildingBlock("c1ccccc1"),
            rmsd=0.0,
            kabsch_rmsd=0.0,
        ),
    ],
)
def different_case_data(request: pytest.FixtureRequest) -> CaseData:
    """A pair of :class:`stk.Molecule` instances and an RMSD."""
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=stk.BuildingBlock(
                smiles=(
                    "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cncc"
                    "c54)cc3c2c1"
                ),
            ),
            mol2=stk.BuildingBlock(
                smiles=(
                    "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cncc"
                    "c54)cc3c2c1"
                ),
            ).with_canonical_atom_ordering(),
            rmsd=0.0,
            kabsch_rmsd=0.0,
        ),
    ],
)
def ordering_case_data(request: pytest.FixtureRequest) -> CaseData:
    """A pair of :class:`stk.Molecule` instances and an RMSD."""
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_centroid(np.array((4, 0, 0))),
            rmsd=0.0,
            kabsch_rmsd=0.0,
        ),
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_position_matrix(
                np.array(
                    [[0.7520009, 0.0, 0.0], [-0.7520009, 0.0, 0.0]],
                )
            ),
            rmsd=0.0,
            kabsch_rmsd=0.0,
        ),
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_position_matrix(
                np.array(
                    [[1.7520009, 0.0, 0.0], [-1.7520009, 0.0, 0.0]],
                )
            ),
            rmsd=1.0,
            kabsch_rmsd=1.0,
        ),
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=stk.BuildingBlock("NCCN")
            .with_rotation_about_axis(
                1.34,
                np.array((0, 0, 1)),
                np.array((0, 0, 0)),
            )
            .with_displacement(np.array((2, 0, 1))),
            rmsd=1.1483491262483871,
            kabsch_rmsd=1.1309858484314543,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=stk.BuildingBlock("CCCCCC")
            .with_rotation_about_axis(
                0.24,
                np.array((1, 0, 1)),
                np.array((0, 0, 0)),
            )
            .with_displacement(np.array((0, 0, 1))),
            rmsd=0.6074500611659212,
            kabsch_rmsd=0.5943193981905652,
        ),
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=stk.BuildingBlock("NCCCN"),
            rmsd=0.8775932363034411,
            kabsch_rmsd=0.8832914099448816,
        ),
        CaseData(
            mol1=stk.BuildingBlock("NCOCN"),
            mol2=stk.BuildingBlock("NCCN"),
            rmsd=1.2590308595881023,
            kabsch_rmsd=1.2678595995702466,
        ),
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=stk.BuildingBlock("NCOCN"),
            rmsd=1.3886830389567653,
            kabsch_rmsd=1.4164013380180553,
        ),
    ],
)
def aligned_case_data(request: pytest.FixtureRequest) -> CaseData:
    """A pair of :class:`stk.Molecule` instances and an RMSD."""
    return request.param
