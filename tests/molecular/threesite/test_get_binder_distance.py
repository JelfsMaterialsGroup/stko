import numpy as np

import stko
from tests.molecular.threesite.case_data import CaseData


def test_get_binder_distance(case_data: CaseData) -> None:
    threesite_analysis = stko.molecule_analysis.DitopicThreeSiteAnalyser()

    distance = threesite_analysis.get_binder_distance(case_data.building_block)
    assert np.isclose(case_data.binder_distance, distance, rtol=0, atol=1e-3)
