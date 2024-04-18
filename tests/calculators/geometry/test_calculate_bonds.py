import numpy as np
import stko


def test_calculate_bonds(case_data):
    """Test :class:`.GeometryAnalyser.calculate_bonds`.

    Parameters
    ----------
        case_data:
            A test case.

    """
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.calculate_bonds(case_data.molecule)
    print(result)
    for pair in result:
        if pair == ("C", "C"):
            continue
        if "H" in pair:
            continue
        assert len(result[pair]) == len(case_data.bonds[pair])
        assert np.allclose(
            result[pair], case_data.bonds[pair], rtol=0, atol=1e-3
        )
