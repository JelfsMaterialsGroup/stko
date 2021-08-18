import stko


def test_zmatrix(case_data):
    """
    Test :meth:`.Converter.get_zmatrix`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        The test case.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_zmatrix(case_data.molecule, case_data.zmatrix)


def _test_zmatrix(molecule, zmatrix):
    """
    Test :meth:`.Converter.get_zmatrix`.

    Parameters
    ----------
    molecule : :class:`stk.Molecule`
        The molecule to test.

    zmatrix : :class:`str`
        The correct Z-matrix string.

    Returns
    -------
    None : :class:`NoneType`

    """

    result = stko.ZMatrix().get_zmatrix(molecule)
    # Printing is useful for debugging.
    print(result)
    assert result == zmatrix
