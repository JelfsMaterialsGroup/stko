import pytest

import stko

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            smiles="C1=CC=CC=C1",
            num_fgs=0,
            factory=stko.functional_groups.ThreeSiteFactory(
                smarts="[Br][C][O]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            smiles="C1=NC=CN=C1",
            num_fgs=2,
            factory=stko.functional_groups.CNCFactory(),
            name=name,
        ),
        lambda name: CaseData(
            smiles="C1=NC=CN=C1",
            num_fgs=0,
            factory=stko.functional_groups.CNNFactory(),
            name=name,
        ),
        lambda name: CaseData(
            smiles="C1=CC=CN=C1",
            num_fgs=1,
            factory=stko.functional_groups.CNCFactory(),
            name=name,
        ),
        lambda name: CaseData(
            smiles="C1=CN=NC=C1",
            num_fgs=0,
            factory=stko.functional_groups.ThreeSiteFactory(
                smarts="[Br][C][O]"
            ),
            name=name,
        ),
        lambda name: CaseData(
            smiles="C1=CN=NC=C1",
            num_fgs=2,
            factory=stko.functional_groups.CNNFactory(),
            name=name,
        ),
        lambda name: CaseData(
            smiles="C1=CN=NC=C1",
            num_fgs=0,
            factory=stko.functional_groups.CNCFactory(),
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
