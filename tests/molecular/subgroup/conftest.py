import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock("C1=CC=CC=C1"),
            sub_group_data={
                "c6_planarity": [2.7518147481201438e-06],
                "c5n1_planarity": [],
                "x5_planarity": [],
                "c#c_angle": [],
            },
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock("C1N=CC(CCC2CCOC2)N=1"),
            sub_group_data={
                "c6_planarity": [],
                "c5n1_planarity": [],
                "x5_planarity": [1.3688005804646254e-06, 0.932064037529801],
                "c#c_angle": [],
            },
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock("C1=CC=C(C=C1)C#CC2=CN=CC=C2"),
            sub_group_data={
                "c6_planarity": [8.41286151020968e-08],
                "c5n1_planarity": [5.678704369238556e-08],
                "x5_planarity": [],
                "c#c_angle": [179.00063441359868],
            },
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
