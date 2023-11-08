import numpy as np
import stko


def test_get_building_block_centroids(case_data):
    """
    Test :class:`.ConstructedAnalyser.get_building_block_centroids`.

    Parameters:

        case_data:
            A test case.

    """

    analyser = stko.molecule_analysis.ConstructedAnalyser()
    result = analyser.get_building_block_centroids(
        case_data.constructed_molecule
    )
    print(result)
    assert len(result) == len(case_data.centroids)
    for i in result:
        print(result[i], case_data.centroids[i])
        assert np.allclose(
            result[i], case_data.centroids[i], rtol=0, atol=1e-3
        )
