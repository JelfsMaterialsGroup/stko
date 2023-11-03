import numpy as np
import stko


def test_get_building_block_centroids(case_data):
    """
    Test :class:`.ConstructedAnalyser.get_building_block_centroids`.

    Parameters:

        case_data:
            A test case.

    """

    assert False

    threesite_analysis = stko.molecule_analysis.DitopicThreeSiteAnalyser()

    distance = threesite_analysis.get_binder_distance(case_data.building_block)
    print(distance)
    assert np.isclose(case_data.binder_distance, distance, rtol=0, atol=1e-3)
