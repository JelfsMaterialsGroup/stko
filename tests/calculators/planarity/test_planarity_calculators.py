import numpy as np

import stko

from .conftest import CaseData


def test_planarity_calculation(case_data: CaseData) -> None:
    calc = stko.PlanarityCalculator()
    test = calc.get_results(
        mol=case_data.molecule,
        plane_atom_ids=case_data.plane_ids,
        deviation_atom_ids=case_data.deviation_ids,
    )
    assert np.isclose(
        test.get_plane_deviation(),
        case_data.plane_deviation,
        atol=1e-4,
    )
    assert np.isclose(
        test.get_plane_deviation_span(),
        case_data.plane_span,
        atol=1e-4,
    )
    assert np.isclose(
        test.get_planarity_parameter(),
        case_data.planarity_parameter,
        atol=1e-4,
    )
