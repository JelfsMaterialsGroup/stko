import numpy as np
import stko


def test_calculate_angles(case_data):
    """
    Test :class:`.GeometryAnalyser.calculate_angles`.

    Parameters:

        case_data:
            A test case.

    """

    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.calculate_angles(case_data.molecule)
    print(result)
    for triple in result:
        if triple == ("C", "C", "C"):
            continue
        if "H" in triple:
            continue
        assert len(result[triple]) == len(case_data.angles[triple])
        assert np.allclose(
            result[triple], case_data.angles[triple], rtol=0, atol=1e-3
        )
