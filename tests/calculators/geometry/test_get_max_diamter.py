import numpy as np

import stko

from .case_data import CaseData


def test_get_max_diameter(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_max_diameter(case_data.molecule)
    assert np.isclose(result, case_data.max_diameter, atol=1e-3, rtol=0)
