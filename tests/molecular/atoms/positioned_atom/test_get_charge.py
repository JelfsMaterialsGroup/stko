from .case_data import CaseData


def test_get_charge(case_data: CaseData) -> None:
    assert case_data.atom.get_charge() == case_data.charge
