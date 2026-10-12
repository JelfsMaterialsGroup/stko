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
            unoptimised_energy=17.81566151445278,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)cc3c2c1"
            ),
            unoptimised_energy=137.86271760336925,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            unoptimised_energy=19.109681124295133,
        ),
        CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            unoptimised_energy=12.140273163773355,
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
            unoptimised_energy=3625.041945631242,
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
            unoptimised_energy=25.785201586676408,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)cc3c2c1"
            ),
            unoptimised_energy=117.7534792935808,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            unoptimised_energy=7.61060374174415,
        ),
        CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            unoptimised_energy=17.198563088658787,
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
            unoptimised_energy=1004.392045956138,
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
            "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)cc3c2c1"
        ),
        stk.BuildingBlock("CCCCCC"),
        stk.BuildingBlock("c1ccccc1"),
    ],
)
def case_etkdg_molecule(request: pytest.FixtureRequest) -> stk.BuildingBlock:
    return request.param
