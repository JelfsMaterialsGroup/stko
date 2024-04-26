from dataclasses import dataclass

import pytest
import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    unoptimised_energy: float


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            unoptimised_energy=18.706050515892986,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)"
                "cc3c2c1"
            ),
            unoptimised_energy=276.0206611549808,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            unoptimised_energy=20.722743438967758,
        ),
        CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            unoptimised_energy=13.516838919531384,
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
                    optimizer=stk.MCHammer(),
                ),
            ),
            unoptimised_energy=5348.367149383393,
        ),
    ],
)
def case_uff_molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            unoptimised_energy=26.518703818643935,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)"
                "cc3c2c1"
            ),
            unoptimised_energy=226.18914087716263,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            unoptimised_energy=7.607569230469989,
        ),
        CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            unoptimised_energy=17.833167064273834,
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
                    optimizer=stk.MCHammer(),
                ),
            ),
            unoptimised_energy=953.310417015842,
        ),
    ],
)
def case_mmff_molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        stk.BuildingBlock("NCCN"),
        stk.BuildingBlock(
            "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)" "cc3c2c1"
        ),
        stk.BuildingBlock("CCCCCC"),
        stk.BuildingBlock("c1ccccc1"),
    ],
)
def case_etkdg_molecule(request: pytest.FixtureRequest) -> stk.BuildingBlock:
    return request.param
