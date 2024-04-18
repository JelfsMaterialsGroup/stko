import pytest
import stk
import stko

_macrocycle = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
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
    ),
)

_square_planar = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.SquarePlanar(
        metals=stk.BuildingBlock(
            smiles="[Pd+2]",
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2)) for i in range(4)
            ),
            position_matrix=[[0, 0, 0]],
        ),
        ligands=stk.BuildingBlock(
            smiles="NBr",
            functional_groups=(stk.PrimaryAminoFactory(),),
        ),
        optimizer=stk.MCHammer(),
    )
)
uff = stko.UFF()
_square_planar_uff = uff.optimize(_square_planar)
_octahedral = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.Octahedral(
        metals=stk.BuildingBlock(
            smiles="[Fe+2]",
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2)) for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        ),
        ligands=stk.BuildingBlock(
            smiles="NBr",
            functional_groups=(stk.PrimaryAminoFactory(),),
        ),
        optimizer=stk.MCHammer(),
    )
)


class CaseData:
    """A test case.

    Attributes:
    ----------
        molecule:
            The molecule to be tested.

        plane_ids:
            The atom ids to define the plane.

        deviation_ids:
            The atom ids to calculate planarity of.

        plane_deviation:
            The plane deviation of the molecule.

        planarity_parameter:
            The planarity parameter of the molecule.

        plane_span:
            The plane deviation span of the molecule.

    """

    def __init__(
        self,
        molecule,
        plane_ids,
        deviation_ids,
        plane_deviation,
        planarity_parameter,
        plane_span,
    ):
        self.molecule = molecule
        self.plane_ids = plane_ids
        self.deviation_ids = deviation_ids
        self.plane_deviation = plane_deviation
        self.planarity_parameter = planarity_parameter
        self.plane_span = plane_span


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=7.011539016627999,
            plane_span=2.881264735947162,
            planarity_parameter=0.7397816598925971,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)"
                "cc3c2c1"
            ),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=34.81901656204046,
            plane_span=3.81412213414559,
            planarity_parameter=0.7871978841119759,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=11.582698264496539,
            plane_span=2.875963748011898,
            planarity_parameter=0.7572794957405952,
        ),
        CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            plane_span=0.0,
            planarity_parameter=0.0,
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
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
                ),
            ),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=15.93962647104854,
            plane_span=2.898937325615976,
            planarity_parameter=0.6223568662643515,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=20.743180343997267,
            plane_span=2.3434639552983,
            planarity_parameter=0.7222254597555785,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6
            ),
            deviation_ids=(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6
            ),
            plane_deviation=1.7532095738358657,
            plane_span=0.5478779918237081,
            planarity_parameter=0.2191511967294832,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6
            ),
            deviation_ids=None,
            plane_deviation=20.694327117503335,
            plane_span=2.3434639552983,
            planarity_parameter=0.7222713549711792,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=None,
            deviation_ids=(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6
            ),
            plane_deviation=1.8020628003297954,
            plane_span=0.5478779918237081,
            planarity_parameter=0.21930239971489357,
        ),
        CaseData(
            molecule=_square_planar,
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            plane_span=0.0,
            planarity_parameter=0.0,
        ),
        CaseData(
            molecule=_square_planar_uff,
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            plane_span=0.0,
            planarity_parameter=0.0,
        ),
        CaseData(
            molecule=_octahedral,
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=15.350980683338454,
            plane_span=5.734353982049391,
            planarity_parameter=1.4635207224067897,
        ),
    ],
)
def case_data(request):
    return request.param
