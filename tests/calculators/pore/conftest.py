import pytest
import stk

from .case_data import CaseData

_complex = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.OctahedralDelta(
        metals=stk.BuildingBlock(
            smiles="[Fe+2]",
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2)) for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
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
            metal_atom_distances={},
            metal_centroid_angles={},
            min_centoid_distance=3.271852034747362,
            avg_centoid_distance=(5.328448699369716, 1.0456647551501264),
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
            metal_atom_distances={
                (0, 40): 14.83956437672709,
                (0, 80): 14.454856249450632,
                (0, 120): 14.837980663967297,
                (40, 80): 14.593027908302588,
                (40, 120): 14.655844535316485,
                (80, 120): 14.599774709094143,
            },
            metal_centroid_angles={
                (0, 40): 110.77513231848364,
                (0, 80): 108.23994750583657,
                (0, 120): 110.55253117949731,
                (40, 80): 109.53287307949009,
                (40, 120): 108.32171727801348,
                (80, 120): 109.40454938926112,
            },
            min_centoid_distance=4.028829338000247,
            avg_centoid_distance=(8.93396048240581, 2.476075586273826),
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
                                for i in range(4)
                            ),
                            position_matrix=[[0, 0, 0]],
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
            metal_atom_distances={
                (0, 1): 10.896944441369676,
                (0, 2): 15.412420828743615,
                (0, 3): 10.898227238903015,
                (0, 4): 10.89758581086942,
                (0, 5): 10.89758581086942,
                (1, 2): 10.898227013424108,
                (1, 3): 15.412420828751983,
                (1, 4): 10.89758569813956,
                (1, 5): 10.89758569813956,
                (2, 3): 10.899510587395719,
                (2, 4): 10.898868771096192,
                (2, 5): 10.898868771096192,
                (3, 4): 10.898868883845243,
                (3, 5): 10.898868883845243,
                (4, 5): 15.412420799323275,
            },
            metal_centroid_angles={
                (0, 1): 89.97683228781653,
                (0, 2): 179.97682458236295,
                (0, 3): 89.99999891399814,
                (0, 4): 89.98841325868929,
                (0, 5): 89.98841325868929,
                (1, 2): 89.99999229454916,
                (1, 3): 179.97683120183706,
                (1, 4): 89.988409950303,
                (1, 5): 89.988409950303,
                (2, 3): 90.02317650363617,
                (2, 4): 90.01158205422921,
                (2, 5): 90.01158205422921,
                (3, 4): 90.0115853652926,
                (3, 5): 90.0115853652926,
                (4, 5): 179.96722969198316,
            },
            min_centoid_distance=6.383408826308891,
            avg_centoid_distance=(8.60220669009043, 1.4865651154918351),
            name=name,
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
