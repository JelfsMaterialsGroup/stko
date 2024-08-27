import numpy as np

import stko
from tests.molecular.threesite.case_data import CaseData


def test_get_binder_angles(case_data: CaseData) -> None:
    threesite_analysis = stko.molecule_analysis.DitopicThreeSiteAnalyser()

    result = threesite_analysis.get_binder_angles(case_data.building_block)
    assert np.isclose(case_data.binder_angles[0], result[0], rtol=0, atol=1e-3)
    assert np.isclose(case_data.binder_angles[1], result[1], rtol=0, atol=1e-3)

    result = threesite_analysis.get_halfbite_angles(case_data.building_block)
    assert np.isclose(case_data.bite_angles[0], result[0], rtol=0, atol=1e-3)
    assert np.isclose(case_data.bite_angles[1], result[1], rtol=0, atol=1e-3)
