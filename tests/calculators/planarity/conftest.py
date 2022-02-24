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

    position_matrix: np.ndarray

    def __init__(self, molecule, unoptimised_energy):

        self.molecule = molecule
        self.unoptimised_energy = unoptimised_energy


@pytest.fixture(
    scope='session',
    params=[
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            unoptimised_energy=18.706050515892986,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)'
                'cc3c2c1'
            ),
            unoptimised_energy=276.0206611549808,
        ),
        CaseData(
            molecule=stk.BuildingBlock('CCCCCC'),
            unoptimised_energy=20.722743438967758,
        ),
        CaseData(
            molecule=stk.BuildingBlock('c1ccccc1'),
            unoptimised_energy=13.516838919531384,
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
            unoptimised_energy=3533.7741683439153,
        ),
    ],
)
def case_data(request):

    return request.param
