import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            zmatrix=(
                "N\n"
                "C 1 1.44\n"
                "C 2 1.531 54.65\n"
                "N 3 1.42 111.07 1 -60.0\n"
                "H 4 3.353 114.1 2 49.59\n"
                "H 5 1.714 52.52 3 -139.34\n"
                "H 6 2.495 49.74 4 65.03\n"
                "H 7 1.846 87.9 5 -7.22\n"
                "H 8 2.377 85.07 6 -58.27\n"
                "H 9 1.818 92.5 7 -59.1\n"
                "H 10 2.99 69.57 8 119.38\n"
                "H 11 1.7510 51.71 9 -161.42"
            ),
        ),
        CaseData(
            molecule=stk.BuildingBlock("BrC#CBr"),
            zmatrix=(
                "Br\n"
                "C 1 1.9\n"
                "C 2 1.211 0.17\n"
                "Br 3 1.892 179.14 1 -171.59"
            ),
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
