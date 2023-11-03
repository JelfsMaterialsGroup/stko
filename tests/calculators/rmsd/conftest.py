import numpy as np
import pytest
import stk
import stko


class CaseData:
    """
    A test case.

    Attributes:
        mol1:
            The first molecule to be tested.

        mol2:
            The second molecule to be tested.

        rmsd:
            The RMSD of the pair of molecules.

    """

    position_matrix: np.ndarray

    def __init__(self, mol1, mol2, rmsd):
        self.mol1 = mol1
        self.mol2 = mol2
        self.rmsd = rmsd


_optimizer = stko.UFF()
_polymer = stk.ConstructedMolecule(
    topology_graph=stk.polymer.Linear(
        building_blocks=(
            stk.BuildingBlock(
                smiles="BrCCBr",
                functional_groups=[stk.BromoFactory()],
            ),
            stk.BuildingBlock(
                smiles="BrCNCCBr",
                functional_groups=[stk.BromoFactory()],
            ),
        ),
        repeating_unit="AB",
        num_repeating_units=2,
        optimizer=stk.MCHammer(),
    ),
)


_cc_molecule = stk.BuildingBlock("[C][C]")


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_centroid(np.array((4, 0, 0))),
            rmsd=0.0,
        ),
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_position_matrix(
                np.array(
                    [[0.7520009, 0.0, 0.0], [-0.7520009, 0.0, 0.0]],
                )
            ),
            rmsd=0.0,
        ),
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_position_matrix(
                np.array(
                    [[1.7520009, 0.0, 0.0], [-1.7520009, 0.0, 0.0]],
                )
            ),
            rmsd=1.0,
        ),
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=_optimizer.optimize(stk.BuildingBlock("NCCN")),
            rmsd=0.24492870054279647,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=stk.BuildingBlock("CCCCCC"),
            rmsd=0.0,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=_optimizer.optimize(stk.BuildingBlock("CCCCCC")),
            rmsd=0.35636491354918015,
        ),
        CaseData(
            mol1=stk.BuildingBlock("c1ccccc1"),
            mol2=_optimizer.optimize(stk.BuildingBlock("c1ccccc1")),
            rmsd=0.02936762392637932,
        ),
        CaseData(
            mol1=_polymer,
            mol2=_optimizer.optimize(_polymer),
            rmsd=2.257184652840373,
        ),
    ],
)
def case_data(request):
    """
    A pair of :class:`stk.Molecule` instances and an RMSD.

    """

    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=_optimizer.optimize(stk.BuildingBlock("NCCN")),
            rmsd=0.20811702035676308,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=_optimizer.optimize(stk.BuildingBlock("CCCCCC")),
            rmsd=0.22563756374632568,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=stk.BuildingBlock("CCCCCC"),
            rmsd=0.0,
        ),
        CaseData(
            mol1=stk.BuildingBlock("c1ccccc1"),
            mol2=_optimizer.optimize(stk.BuildingBlock("c1ccccc1")),
            rmsd=0.029156836455717483,
        ),
        CaseData(
            mol1=_polymer,
            mol2=_optimizer.optimize(_polymer),
            rmsd=2.1395034910834867,
        ),
    ],
)
def ignore_h_case_data(request):
    """
    A pair of :class:`stk.Molecule` instances and an RMSD.

    """

    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=stk.BuildingBlock("CCCCCC"),
            rmsd=0.0,
        ),
        CaseData(
            mol1=stk.BuildingBlock("CCCCCC"),
            mol2=stk.BuildingBlock("c1ccccc1"),
            rmsd=0.0,
        ),
    ],
)
def different_case_data(request):
    """
    A pair of :class:`stk.Molecule` instances and an RMSD.

    """

    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=_polymer,
            mol2=_polymer.with_canonical_atom_ordering(),
            rmsd=0.0,
        ),
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
        ),
    ],
)
def ordering_case_data(request):
    """
    A pair of :class:`stk.Molecule` instances and an RMSD.

    """

    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_centroid(np.array((4, 0, 0))),
            rmsd=0.0,
        ),
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_position_matrix(
                np.array(
                    [[0.7520009, 0.0, 0.0], [-0.7520009, 0.0, 0.0]],
                )
            ),
            rmsd=0.0,
        ),
        CaseData(
            mol1=_cc_molecule,
            mol2=_cc_molecule.with_position_matrix(
                np.array(
                    [[1.7520009, 0.0, 0.0], [-1.7520009, 0.0, 0.0]],
                )
            ),
            rmsd=1.0,
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
            rmsd=1.1309858484314543,
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
            rmsd=0.5943193981905652,
        ),
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=stk.BuildingBlock("NCCCN"),
            rmsd=0.8832914099448816,
        ),
        CaseData(
            mol1=stk.BuildingBlock("NCOCN"),
            mol2=stk.BuildingBlock("NCCN"),
            rmsd=1.2678595995702466,
        ),
        CaseData(
            mol1=stk.BuildingBlock("NCCN"),
            mol2=stk.BuildingBlock("NCOCN"),
            rmsd=1.3921770318522637,
        ),
    ],
)
def aligned_case_data(request):
    """
    A pair of :class:`stk.Molecule` instances and an RMSD.

    """

    return request.param
