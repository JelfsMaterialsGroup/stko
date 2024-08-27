import numpy as np

import stko
from tests.molecular.threesite.case_data import CaseData


def test_get_adjacent_centroids(case_data: CaseData) -> None:
    threesite_analysis = stko.molecule_analysis.DitopicThreeSiteAnalyser()

    centroids = threesite_analysis.get_adjacent_centroids(
        case_data.building_block
    )
    assert np.linalg.norm(centroids[0] - centroids[1]) > 0.1  # noqa: PLR2004
