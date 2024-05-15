import stk

from tests.molecular.functional_group.case_data import CaseData


def test_threesite(case_data: CaseData) -> None:
    num_fgs = stk.BuildingBlock(
        smiles=case_data.smiles,
        functional_groups=case_data.factory,
    ).get_num_functional_groups()

    assert case_data.num_fgs == num_fgs
