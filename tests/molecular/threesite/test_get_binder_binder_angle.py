import numpy as np

import stko
from tests.molecular.threesite.case_data import CaseData


def test_get_binder_binder_angle(case_data: CaseData) -> None:
    threesite_analysis = stko.molecule_analysis.DitopicThreeSiteAnalyser()

    result = threesite_analysis.get_binder_binder_angle(
        case_data.building_block
    )
    assert np.isclose(
        case_data.binder_binder_angle,
        result,
        rtol=0,
        atol=1e-3,
    )
