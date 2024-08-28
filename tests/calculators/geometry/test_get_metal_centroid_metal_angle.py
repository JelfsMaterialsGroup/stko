import numpy as np

import stko

from .case_data import CaseData


def test_get_metal_centroid_metal_angle(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_metal_centroid_metal_angle(
        case_data.molecule,
        metal_atom_nos=(26, 46),
    )
    assert len(result) == len(case_data.metal_centroid_angles)
    for i in result:
        assert np.isclose(
            result[i], case_data.metal_centroid_angles[i], atol=1e-3, rtol=0
        )
