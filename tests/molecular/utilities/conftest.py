from dataclasses import dataclass

import pytest
import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    molecules: stk.Molecule

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
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
