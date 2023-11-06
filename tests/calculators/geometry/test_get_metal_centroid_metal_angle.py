import numpy as np
import stko


def test_get_metal_centroid_metal_angle(case_data):
    """
    Test :class:`.GeometryAnalyser.get_metal_centroid_metal_angle`.

    Parameters:

        case_data:
            A test case.

    """

    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_metal_centroid_metal_angle(
        case_data.molecule,
        metal_atom_nos=(26, 46),
    )
    print(result)
    assert len(result) == len(case_data.metal_centroid_angles)
    for i in result:
        print(i, result[i])
        assert np.isclose(
            result[i], case_data.metal_centroid_angles[i], atol=1e-3, rtol=0
        )
