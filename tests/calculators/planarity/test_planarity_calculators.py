import json
import pathlib

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
    test_dict = {
        "atom_plane": [float(i) for i in test.get_atom_plane()],
        "planarity_parameter": test.get_planarity_parameter(),
        "plane_deviation": test.get_plane_deviation(),
        "plane_deviation_span": test.get_plane_deviation_span(),
    }

    parent = pathlib.Path(__file__).resolve().parent
    saved_dir = parent / "planarity_jsons"
    saved_json = saved_dir / f"{case_data.name}.json"
    if not saved_json.exists():
        with saved_json.open("w") as f:
            json.dump(test_dict, f, indent=4)

    with saved_json.open("r") as f:
        known = json.load(f)

    assert len(test_dict) == len(known)

    for i, test_value in test_dict.items():
        assert np.allclose(test_value, known[str(i)], rtol=0, atol=1e-3)
