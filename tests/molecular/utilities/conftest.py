from dataclasses import dataclass

import pytest
import stk

import stko


@dataclass(frozen=True, slots=True)
class CaseData:
    molecules: list[stk.Molecule]
    name: str


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecules=[stk.BuildingBlock("C1=CC=CC=C1")],
            name=name,
        ),
        lambda name: CaseData(
            molecules=[
                stk.BuildingBlock("C1=CC=CC=C1"),
                stk.BuildingBlock("C1N=CC(CCC2CCOC2)N=1"),
            ],
            name=name,
        ),
        lambda name: CaseData(
            molecules=[
                stk.BuildingBlock("C1=CC=CC=C1"),
                stk.BuildingBlock("C1=CC=C(C=C1)C#CC2=CN=CC=C2"),
            ],
            name=name,
        ),
        lambda name: CaseData(
            molecules=[
                stk.BuildingBlock("C1=CC=CC=C1"),
                stk.BuildingBlock("C1N=CC(CCC2CCOC2)N=1"),
                stk.BuildingBlock("C1=CC=C(C=C1)C#CC2=CN=CC=C2"),
            ],
            name=name,
        ),
        lambda name: CaseData(
            molecules=[
                i[0]
                for i in stko.molecular_utilities.separate_molecule(
                    stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(
                            host=stk.BuildingBlock("NCCN"),
                            guests=stk.host_guest.Guest(
                                building_block=stk.BuildingBlock("NCCN"),
                                displacement=(5, 5, 5),
                            ),
                        ),
                    )
                )
            ],
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )


@dataclass(frozen=True, slots=True)
class TogetherCaseData:
    molecule: stk.Molecule
    name: str


@pytest.fixture(
    scope="session",
    params=(
        lambda name: TogetherCaseData(
            molecule=stk.BuildingBlock("C1=CC=CC=C1"),
            name=name,
        ),
        lambda name: TogetherCaseData(
            molecule=stko.molecular_utilities.merge_stk_molecules(
                [
                    stk.BuildingBlock("C1=CC=CC=C1"),
                    stk.BuildingBlock("C1N=CC(CCC2CCOC2)N=1"),
                    stk.BuildingBlock("C1=CC=C(C=C1)C#CC2=CN=CC=C2"),
                ]
            ),
            name=name,
        ),
        lambda name: TogetherCaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock("NCCN"),
                    guests=stk.host_guest.Guest(
                        building_block=stk.BuildingBlock("NCCN"),
                        displacement=(5, 5, 5),
                    ),
                ),
            ),
            name=name,
        ),
        lambda name: TogetherCaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock("C1N=CC(CCC2CCOC2)N=1"),
                    guests=stk.host_guest.Guest(
                        building_block=stk.BuildingBlock("NCCN"),
                        displacement=(5, 5, 5),
                    ),
                ),
            ),
            name=name,
        ),
        lambda name: TogetherCaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock("C1N=CC(CCC2CCOC2)N=1"),
                    guests=(
                        stk.host_guest.Guest(
                            building_block=stk.BuildingBlock("NCCN"),
                            displacement=(5, 5, 5),
                        ),
                        stk.host_guest.Guest(
                            building_block=stk.BuildingBlock("C1=CC=CC=C1"),
                            displacement=(5, 5, 5),
                        ),
                    ),
                ),
            ),
            name=name,
        ),
    ),
)
def together_case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
