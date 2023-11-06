import stko


def test_get_metal_distances(case_data):
    """
    Test :class:`.PoreAnalyser.get_metal_distances`.

    Parameters:

        case_data:
            A test case.

    """

    analyser = stko.molecule_analysis.PoreAnalyser()

    result = analyser.get_metal_distances(
        case_data.molecule,
        metal_atom_nos=(26, 46),
    )
    print(result)
    assert len(result) == len(case_data.metal_atom_distances)
    for i in result:
        print(i, result[i])
        assert result[i] == case_data.metal_atom_distances[i]
