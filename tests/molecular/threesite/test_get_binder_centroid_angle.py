import numpy as np
import stko


def test_get_binder_centroid_angle(case_data):
    """
    Test :class:`.DitopicThreeSiteAnalyser.get_binder_centroid_angle`.

    Parameters:

        case_data:
            A test case.

    """

    threesite_analysis = stko.molecule_analysis.DitopicThreeSiteAnalyser()

    result = threesite_analysis.get_binder_centroid_angle(
        case_data.building_block
    )
    print(result)
    assert np.isclose(
        case_data.binder_centroid_angle,
        result,
        rtol=0,
        atol=1e-3,
    )
