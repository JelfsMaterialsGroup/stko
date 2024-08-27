import numpy as np

import stko

from .case_data import CaseData


def test_calculate_angles(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.calculate_angles(case_data.molecule)
    for triple in result:
        if triple == ("C", "C", "C"):
            continue
        if "H" in triple:
            continue
        assert len(result[triple]) == len(case_data.angles[triple])
        assert np.allclose(
            result[triple], case_data.angles[triple], rtol=0, atol=1e-3
        )
