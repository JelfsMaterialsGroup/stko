import stko


def test_get_building_block_atom_ids(case_data):
    """
    Test :class:`.ConstructedAnalyser.get_building_block_atom_ids`.

    Parameters:

        case_data:
            A test case.

    """

    analyser = stko.molecule_analysis.ConstructedAnalyser()
    result = analyser.get_building_block_atom_ids(
        case_data.constructed_molecule
    )
    print(result)
    assert len(result) == len(case_data.atom_ids)
    for i in result:
        print(result[i], case_data.atom_ids[i])
        assert result[i] == case_data.atom_ids[i]
