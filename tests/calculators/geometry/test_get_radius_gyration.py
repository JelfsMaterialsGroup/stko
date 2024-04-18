import numpy as np
import stko


def test_get_radius_gyration(case_data):
    """Test :class:`.GeometryAnalyser.get_radius_gyration`.

    Parameters
    ----------
        case_data:
            A test case.

    """
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_radius_gyration(case_data.molecule)
    print(result)
    assert np.isclose(result, case_data.radius_gyration, atol=1e-3, rtol=0)
