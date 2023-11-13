import numpy as np
import stko


def test_calculate_torsions(case_data):
    """
    Test :class:`.GeometryAnalyser.calculate_torsions`.

    Parameters:

        case_data:
            A test case.

    """

    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.calculate_torsions(case_data.molecule)
    print(result)
    for four in result:
        if four == ("C", "C", "C", "C"):
            continue
        if "H" in four:
            continue
        assert len(result[four]) == len(case_data.torsions[four])
        assert np.allclose(
            result[four], case_data.torsions[four], rtol=0, atol=1e-3
        )
