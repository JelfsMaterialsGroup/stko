import numpy as np

import stko
from tests.molecular.constructed.case_data import CaseData


def test_get_building_block_centroids(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.ConstructedAnalyser()
    result = analyser.get_building_block_centroids(
        case_data.constructed_molecule
    )
    assert len(result) == len(case_data.centroids)
    for i in result:
        assert np.allclose(
            result[i], case_data.centroids[i], rtol=0, atol=1e-3
        )
