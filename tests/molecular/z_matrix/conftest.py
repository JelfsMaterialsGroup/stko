import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            zmatrix=(
                "N\nC 1 1.46\nC 2 1.51 53.82\nN 3 1.412 112.72 1 -60.01\nH 4 3.413 116.31 2 49.41\nH 5 1.744 51.39 3 -134.65\nH 6 2.475 49.76 4 66.66\nH 7 1.86 88.57 5 -7.45\nH 8 2.377 85.03 6 -59.05\nH 9 1.838 93.26 7 -58.76\nH 10 2.949 70.39 8 117.46\nH 11 1.7510 50.5 9 -162.37"
            ),
        ),
        CaseData(
            molecule=stk.BuildingBlock("BrC#CBr"),
            zmatrix=(
                "Br\nC 1 1.89\nC 2 1.231 0.28\nBr 3 1.892 179.03 1 -161.48"
            ),
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
