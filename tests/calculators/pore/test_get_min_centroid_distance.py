import numpy as np
import stko


def test_get_min_centroid_distance(case_data):
    """
    Test :class:`.PoreAnalyser.get_min_centroid_distance`.

    Parameters:

        case_data:
            A test case.

    """

    analyser = stko.molecule_analysis.PoreAnalyser()

    result = analyser.get_min_centroid_distance(case_data.molecule)
    assert np.isclose(
        result, case_data.min_centoid_distance, atol=1e-3, rtol=0
    )
