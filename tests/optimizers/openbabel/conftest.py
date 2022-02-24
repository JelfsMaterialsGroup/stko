import pytest
import stk


class CaseData:
    """
    A test case.

    Attributes:
        molecule:
            The molecule to be tested.

        unoptimised_energy:
            The energy of the molecule from stk generation.

    """

    def __init__(self, molecule, unoptimised_energy):

        self.molecule = molecule
        self.unoptimised_energy = unoptimised_energy


@pytest.fixture(
    scope='session',
    params=[
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            unoptimised_energy=98.68403256936926,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)'
                'cc3c2c1'
            ),
            unoptimised_energy=2032.8111567743895,
        ),
        CaseData(
            molecule=stk.BuildingBlock('CCCCCC'),
            unoptimised_energy=141.44622279628743,
        ),
        CaseData(
            molecule=stk.BuildingBlock('c1ccccc1'),
            unoptimised_energy=56.73128282588534,
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
                    optimizer=stk.MCHammer(),
                ),
            ),
            unoptimised_energy=15002.946293854524,
        ),
    ],
)
def case_molecule(request):
    return request.param
