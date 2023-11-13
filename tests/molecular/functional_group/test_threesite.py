import stk


def test_threesite(case_data):
    """
    Test :class:`.ThreeSiteFactory`.

    Parameters:

        case_data:
            A test case.

    """

    num_fgs = stk.BuildingBlock(
        smiles=case_data.smiles,
        functional_groups=case_data.factory,
    ).get_num_functional_groups()

    assert case_data.num_fgs == num_fgs
