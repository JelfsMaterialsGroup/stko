import pytest
import stk


class CaseData:
    """
    A test case.

    Attributes:
        molecule:
            The molecule to be tested.

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
        plane_deviation,
        planarity_parameter,
        plane_deviation_span,
    ):
        self.molecule = molecule
        self.plane_deviation = plane_deviation
        self.planarity_parameter = planarity_parameter
        self.plane_deviation_span = plane_deviation_span


@pytest.fixture(
    scope='session',
    params=[
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            plane_deviation = 0.0,
            planarity_parameter = 0.0,
            plane_deviation_span = 0.0,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)'
                'cc3c2c1'
            ),
            plane_deviation = 0.0,
            planarity_parameter = 0.0,
            plane_deviation_span = 0.0,
        ),
        CaseData(
            molecule=stk.BuildingBlock('CCCCCC'),
            plane_deviation = 0.0,
            planarity_parameter = 0.0,
            plane_deviation_span = 0.0,
        ),
        CaseData(
            molecule=stk.BuildingBlock('c1ccccc1'),
            plane_deviation = 0.0,
            planarity_parameter = 0.0,
            plane_deviation_span = 0.0,
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
            plane_deviation = 0.0,
            planarity_parameter = 0.0,
            plane_deviation_span = 0.0,
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
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
            ),
            plane_deviation = 0.0,
            planarity_parameter = 0.0,
            plane_deviation_span = 0.0,
        ),
    ],
)
def case_data(request):

    return request.param
