import pytest
import stk
import stko

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            constructed_molecule=stk.BuildingBlock(
                smiles="C1=CN=CC=N1",
                functional_groups=stko.functional_groups.ThreeSiteFactory(
                    smarts="[#6]~[#7X2]~[#6]"
                ),
            ),
            name=name,
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
