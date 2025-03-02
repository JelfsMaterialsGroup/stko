from dataclasses import dataclass

import pytest
import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    opt_energy: float


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            molecule=stk.BuildingBlock("C"),
            opt_energy=0.22629718073104876,
        ),
    ],
)
def case_molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param
