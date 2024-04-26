from .case_data import CaseData


def test_get_atomic_number(case_data: CaseData) -> None:
    assert case_data.atom.get_atomic_number() == 0
