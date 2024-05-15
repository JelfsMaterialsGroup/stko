from .case_data import CaseData


def test_get_id(case_data: CaseData) -> None:
    assert case_data.atom.get_id() == case_data.id
