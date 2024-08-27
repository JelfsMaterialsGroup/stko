import stko
from tests.molecular.constructed.case_data import CaseData


def test_get_building_block_atom_ids(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.ConstructedAnalyser()
    result = analyser.get_building_block_atom_ids(
        case_data.constructed_molecule
    )
    assert len(result) == len(case_data.atom_ids)
    for i in result:
        assert result[i] == case_data.atom_ids[i]
