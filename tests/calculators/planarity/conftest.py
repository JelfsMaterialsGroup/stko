import pytest
import stk


_macrocycle = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
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
    ),
)


class CaseData:
    """
    A test case.

    Attributes:
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

        plane_deviation_span:
            The plane deviation span of the molecule.

    """

    def __init__(
        self,
        molecule,
        plane_ids,
        deviation_ids,
        plane_deviation,
        planarity_parameter,
        plane_deviation_span,
    ):
        self.molecule = molecule
        self.plane_ids = plane_ids
        self.deviation_ids = deviation_ids
        self.plane_deviation = plane_deviation
        self.planarity_parameter = planarity_parameter
        self.plane_deviation_span = plane_deviation_span


@pytest.fixture(
    scope='session',
    params=[
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            planarity_parameter=0.0,
            plane_deviation_span=0.0,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)'
                'cc3c2c1'
            ),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            planarity_parameter=0.0,
            plane_deviation_span=0.0,
        ),
        CaseData(
            molecule=stk.BuildingBlock('CCCCCC'),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            planarity_parameter=0.0,
            plane_deviation_span=0.0,
        ),
        CaseData(
            molecule=stk.BuildingBlock('c1ccccc1'),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            planarity_parameter=0.0,
            plane_deviation_span=0.0,
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
                ),
            ),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            planarity_parameter=0.0,
            plane_deviation_span=0.0,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            planarity_parameter=0.0,
            plane_deviation_span=0.0,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=(
                i.get_id() for i in _macrocycle.get_atoms()
                if i.get_atomic_num() == 6
            ),
            deviation_ids=(
                i.get_id() for i in _macrocycle.get_atoms()
                if i.get_atomic_num() == 6
            ),
            plane_deviation=0.0,
            planarity_parameter=0.0,
            plane_deviation_span=0.0,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=(
                i.get_id() for i in _macrocycle.get_atoms()
                if i.get_atomic_num() == 6
            ),
            deviation_ids=None,
            plane_deviation=0.0,
            planarity_parameter=0.0,
            plane_deviation_span=0.0,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=None,
            deviation_ids=(
                i.get_id() for i in _macrocycle.get_atoms()
                if i.get_atomic_num() == 6
            ),
            plane_deviation=0.0,
            planarity_parameter=0.0,
            plane_deviation_span=0.0,
        ),
    ],
)
def case_data(request):

    return request.param
