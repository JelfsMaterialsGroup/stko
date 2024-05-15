import numpy as np
import pytest
import stk

from .case_data import CaseData

_complex = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.OctahedralDelta(
        metals=stk.BuildingBlock(
            smiles="[Fe+2]",
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2)) for _ in range(6)
            ),
            position_matrix=np.array([[0, 0, 0]]),
        ),
        ligands=stk.BuildingBlock(
            smiles="C1=NC(C=NBr)=CC=C1",
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#35]",
                    bonders=(1,),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]",
                    bonders=(1,),
                    deleters=(),
                ),
            ],
        ),
        optimizer=stk.MCHammer(),
    ),
)


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="NCCN",
                            functional_groups=[stk.PrimaryAminoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles="O=CC(C=O)C=O",
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                    ),
                    optimizer=stk.MCHammer(),
                ),
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M4L6TetrahedronSpacer(
                    building_blocks={
                        stk.BuildingBlock.init_from_molecule(
                            molecule=_complex,
                            functional_groups=[stk.BromoFactory()],
                        ): (0, 1, 2, 3),
                        stk.BuildingBlock(
                            smiles=("C1=CC(=CC=C1C2=CC=C(C=C2)Br)Br"),
                            functional_groups=[stk.BromoFactory()],
                        ): (4, 5, 6, 7, 8, 9),
                    },
                    optimizer=stk.MCHammer(),
                ),
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M6L12Cube(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="[Pd+2]",
                            functional_groups=(
                                stk.SingleAtom(stk.Pd(0, charge=2))
                                for _ in range(4)
                            ),
                            position_matrix=np.array([[0, 0, 0]]),
                        ),
                        stk.BuildingBlock(
                            smiles=(
                                "C1=NC=CC(C2=CC=CC(C3=C" "C=NC=C3)=C2)=C1"
                            ),
                            functional_groups=[
                                stk.SmartsFunctionalGroupFactory(
                                    smarts="[#6]~[#7X2]~[#6]",
                                    bonders=(1,),
                                    deleters=(),
                                ),
                            ],
                        ),
                    ),
                    optimizer=stk.Collapser(),
                ),
            ),
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
