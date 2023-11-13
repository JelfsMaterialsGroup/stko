import numpy as np
import stko


def test_get_adjacent_centroids(case_data):
    """
    Test :class:`.DitopicThreeSiteAnalyser.get_adjacent_centroids`.

    Parameters:

        case_data:
            A test case.

    """

    threesite_analysis = stko.molecule_analysis.DitopicThreeSiteAnalyser()

    centroids = threesite_analysis.get_adjacent_centroids(
        case_data.building_block
    )
    print(centroids, np.linalg.norm(centroids[0] - centroids[1]))
    assert np.linalg.norm(centroids[0] - centroids[1]) > 0.1
