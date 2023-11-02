import pytest
import stk
import stko

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            building_block=stk.BuildingBlock(
                smiles="C1=CN=CC=N1",
                functional_groups=stko.functional_groups.ThreeSiteFactory(
                    smarts="[#6]~[#7X2]~[#6]"
                ),
            ),
            binder_distance=2.717,
            name=name,
        ),
        lambda name: CaseData(
            building_block=stk.BuildingBlock(
                smiles="C1=CC(=CN=C1)C2=CC=C(C=C2)C3=CN=CC=C3",
                functional_groups=stko.functional_groups.ThreeSiteFactory(
                    smarts="[#6]~[#7X2]~[#6]"
                ),
            ),
            binder_distance=10.105,
            name=name,
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
