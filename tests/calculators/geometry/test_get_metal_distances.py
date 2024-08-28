import numpy as np

import stko

from .case_data import CaseData


def test_get_metal_distances(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_metal_distances(
        case_data.molecule,
        metal_atom_nos=(26, 46),
    )
    assert len(result) == len(case_data.metal_atom_distances)
    for i in result:
        assert np.isclose(
            result[i], case_data.metal_atom_distances[i], rtol=0, atol=1e-2
        )
