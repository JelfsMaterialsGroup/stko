def test_get_atomic_number(case_data):
    """
    Test :meth:`.Du.get_atomic_number`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case, containing the atom to test and its correct atomic
        number.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.atom.get_atomic_number() == 0
