import numpy as np
import pytest
import stk

from .case_data import CaseData

_reaction_factory = stk.DativeReactionFactory(
    stk.GenericReactionFactory(
        bond_orders={
            frozenset(
                {
                    stk.GenericFunctionalGroup,
                    stk.SingleAtom,
                }
            ): 9,
        },
    ),
)

_fe_complex = stk.ConstructedMolecule(
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
    ),
)


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            cage=stk.ConstructedMolecule(
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
                ),
            ),
            num_ligands=1,
            metal_atom_nos=(46,),
            bb_smiles=(
                r"C1=N\CC/N=C\[C@H]2/C=N\CC/N=C/[C@H]3/C=N\CC/N=C"
                r"\[C@@H]/1/C=N/CC/N=C\[C@@H](/C=N\CC/N=C/2)/C=N/CC/N=C\3",
            ),
            name=name,
        ),
        lambda name: CaseData(
            cage=stk.ConstructedMolecule(
                topology_graph=stk.cage.M2L4Lantern(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="[Pd+2]",
                            functional_groups=(
                                stk.SingleAtom(stk.Pd(0, charge=2))
                                for _ in range(4)
                            ),
                            position_matrix=np.array([[0.0, 0.0, 0.0]]),
                        ),
                        stk.BuildingBlock(
                            smiles="C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1",
                            functional_groups=[
                                stk.SmartsFunctionalGroupFactory(
                                    smarts="[#6]~[#7X2]~[#6]",
                                    bonders=(1,),
                                    deleters=(),
                                ),
                            ],
                        ),
                    ),
                    reaction_factory=_reaction_factory,
                ),
            ),
            num_ligands=4,
            metal_atom_nos=(46,),
            bb_smiles=("c1cc(-c2ccncc2)cc(-c2ccncc2)c1",),
            name=name,
        ),
        lambda name: CaseData(
            cage=stk.ConstructedMolecule(
                topology_graph=stk.cage.M6L12Cube(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="[Pd+2]",
                            functional_groups=(
                                stk.SingleAtom(stk.Pd(0, charge=2))
                                for i in range(4)
                            ),
                            position_matrix=np.array([[0.0, 0.0, 0.0]]),
                        ),
                        stk.BuildingBlock(
                            smiles=("C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1"),
                            functional_groups=[
                                stk.SmartsFunctionalGroupFactory(
                                    smarts="[#6]~[#7X2]~[#6]",
                                    bonders=(1,),
                                    deleters=(),
                                ),
                            ],
                        ),
                    ),
                    reaction_factory=_reaction_factory,
                ),
            ),
            num_ligands=12,
            metal_atom_nos=(46,),
            bb_smiles=("c1cc(-c2ccncc2)cc(-c2ccncc2)c1",),
            name=name,
        ),
        lambda name: CaseData(
            cage=stk.ConstructedMolecule(
                topology_graph=stk.cage.M4L4Tetrahedron(
                    building_blocks={
                        stk.BuildingBlock.init_from_molecule(
                            molecule=_fe_complex,
                            functional_groups=[stk.BromoFactory()],
                        ): (0, 1, 2, 3),
                        stk.BuildingBlock(
                            smiles=("C1=C(C=C(C=C1Br)Br)Br"),
                            functional_groups=[stk.BromoFactory()],
                        ): (4, 5, 6, 7),
                    }
                ),
            ),
            num_ligands=4,
            metal_atom_nos=(26,),
            bb_smiles=(
                r"C(=N/c1cc(/N=C/c2ccccn2)cc(/N=C/c2ccccn2)c1)\c1ccccn1",
            ),
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
