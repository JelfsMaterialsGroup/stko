import numpy as np
import stko


def test_get_min_atom_atom_distance(case_data):
    """Test :class:`.GeometryAnalyser.get_min_atom_atom_distance`.

    Parameters
    ----------
        case_data:
            A test case.

    """
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_min_atom_atom_distance(case_data.molecule)
    print(result)
    assert np.isclose(
        result, case_data.min_atom_atom_distance, atol=1e-3, rtol=0
    )
