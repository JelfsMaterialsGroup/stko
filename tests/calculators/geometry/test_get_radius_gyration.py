import numpy as np

import stko

from .case_data import CaseData


def test_get_radius_gyration(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_radius_gyration(case_data.molecule)
    assert np.isclose(result, case_data.radius_gyration, atol=1e-3, rtol=0)
