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
            binder_centroid_angle=179.55163,
            binder_angles=(179.82673252797306, 178.99117497251564),
            binder_adjacent_torsion=179.98930945981346,
            binder_binder_angle=178.81790750306536,
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
            binder_centroid_angle=179.46618,
            binder_angles=(130.57317931302708, 132.5699156645761),
            binder_adjacent_torsion=-179.7767900679974,
            binder_binder_angle=177.99497293128238,
            name=name,
        ),
        lambda name: CaseData(
            building_block=stk.BuildingBlock(
                smiles="C1=CC(=CC(=C1)C2=CC=NC=C2)C3=CC=NC=C3",
                functional_groups=stko.functional_groups.ThreeSiteFactory(
                    smarts="[#6]~[#7X2]~[#6]"
                ),
            ),
            binder_distance=9.8893,
            binder_centroid_angle=149.0484,
            binder_angles=(152.4051397865058, 153.24145338156575),
            binder_adjacent_torsion=1.1263317393587988,
            binder_binder_angle=125.64880123145808,
            name=name,
        ),
        lambda name: CaseData(
            building_block=stk.BuildingBlock(
                smiles="C1C=C(C2=CC3=C(OC4=C3C=C(C3C=CN=CC=3)C=C4)C=C2)C=CN=1",
                functional_groups=stko.functional_groups.ThreeSiteFactory(
                    smarts="[#6]~[#7X2]~[#6]"
                ),
            ),
            binder_distance=11.778023252068746,
            binder_centroid_angle=132.7736587705255,
            binder_angles=(130.41599653682567, 135.4268203184301),
            binder_adjacent_torsion=-2.7155067603901997,
            binder_binder_angle=85.87201566260612,
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
