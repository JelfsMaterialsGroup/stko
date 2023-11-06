import numpy as np
import stko


def test_get_avg_centroid_distance(case_data):
    """
    Test :class:`.PoreAnalyser.get_avg_centroid_distance`.

    Parameters:

        case_data:
            A test case.

    """

    analyser = stko.molecule_analysis.PoreAnalyser()

    result = analyser.get_avg_centroid_distance(case_data.molecule)
    print(result)
    assert np.isclose(
        result[0], case_data.avg_centoid_distance[0], atol=1e-3, rtol=0
    )
    assert np.isclose(
        result[1], case_data.avg_centoid_distance[1], atol=1e-3, rtol=0
    )
