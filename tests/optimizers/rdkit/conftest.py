import pytest
import numpy as np
import stk


class CaseData:
    """
    A test case.

    Attributes:
        molecule:
            The molecule to be tested.

        unoptimised_energy:
            The energy of the molecule from stk generation.

        optimised_energy:
            The energy of the molecule after optimisation.

    """

    position_matrix: np.ndarray

    def __init__(
        self,
        molecule,
        unoptimised_energy,
        optimised_energy,
    ):

        self.molecule = molecule
        self.unoptimised_energy = unoptimised_energy
        self.optimised_energy = optimised_energy


@pytest.fixture(
    scope='session',
    params=[
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            unoptimised_energy=18.706050515892986,
            optimised_energy=4.873693301322547,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)'
                'cc3c2c1'
            ),
            unoptimised_energy=276.0206611549808,
            optimised_energy=111.21858852673425,
        ),
        CaseData(
            molecule=stk.BuildingBlock('CCCCCC'),
            unoptimised_energy=20.722743438967758,
            optimised_energy=4.424447413069915,
        ),
        CaseData(
            molecule=stk.BuildingBlock('c1ccccc1'),
            unoptimised_energy=13.516838919531384,
            optimised_energy=10.544731800632466,
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
            optimised_energy=15.634852399099843,
        ),
    ],
)
def case_uff_molecule(request):
    """
    A :class:`.Molecule` instance.

    """

    return request.param


@pytest.fixture(
    scope='session',
    params=[
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            unoptimised_energy=26.518703818643935,
            optimised_energy=15.271794627980917,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)'
                'cc3c2c1'
            ),
            unoptimised_energy=226.18914087716263,
            optimised_energy=92.5337600166967,
        ),
        CaseData(
            molecule=stk.BuildingBlock('CCCCCC'),
            unoptimised_energy=7.607569230469989,
            optimised_energy=-4.647071723643174,
        ),
        CaseData(
            molecule=stk.BuildingBlock('c1ccccc1'),
            unoptimised_energy=17.833167064273834,
            optimised_energy=16.22696681429067,
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
            unoptimised_energy=977.3866956352667,
            optimised_energy=8.401610513106467,
        ),
    ],
)
def case_mmff_molecule(request):
    """
    A :class:`.Molecule` instance.

    """

    return request.param


@pytest.fixture(
    scope='session',
    params=[
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock(
            'C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)'
            'cc3c2c1'
        ),
        stk.BuildingBlock('CCCCCC'),
        stk.BuildingBlock('c1ccccc1'),
    ],
)
def case_etkdg_molecule(request):
    """
    A :class:`.Molecule` instance.

    """

    return request.param
