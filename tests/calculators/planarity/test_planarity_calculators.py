import stko
import stk


def test_plane_deviation(case_data):
    calc = stko.PlanarityCalculator()
    result = calc.get_results(
        mol=case_data.molecule,
        plane_atom_ids=case_data.plane_ids,
        deviation_atom_ids=case_data.deviation_ids,
    )
    assert False


def test_plane_deviation_span(case_data):
    calc = stko.PlanarityCalculator()
    result = calc.get_results(
        mol=case_data.molecule,
        plane_atom_ids=case_data.plane_ids,
        deviation_atom_ids=case_data.deviation_ids,
    )
    assert False


def test_planarity_parameter(case_data):
    calc = stko.PlanarityCalculator()
    result = calc.get_results(
        mol=case_data.molecule,
        plane_atom_ids=case_data.plane_ids,
        deviation_atom_ids=case_data.deviation_ids,
    )
    assert False
