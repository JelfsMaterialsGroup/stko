import numpy as np
import stko


def test_get_max_diameter(case_data):
    """
    Test :class:`.GeometryAnalyser.get_max_diameter`.

    Parameters:

        case_data:
            A test case.

    """

    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_max_diameter(case_data.molecule)
    print(result)
    assert np.isclose(result, case_data.max_diameter, atol=1e-3, rtol=0)
