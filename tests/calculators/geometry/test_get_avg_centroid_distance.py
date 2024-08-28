import numpy as np

import stko

from .case_data import CaseData


def test_get_avg_centroid_distance(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_avg_centroid_distance(case_data.molecule)
    assert np.isclose(
        result[0], case_data.avg_centoid_distance[0], atol=1e-3, rtol=0
    )
    assert np.isclose(
        result[1], case_data.avg_centoid_distance[1], atol=1e-3, rtol=0
    )
