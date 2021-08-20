def test_get_position(case_data):
    """
    Test :meth:`.PositionedAtom.get_position`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case, containing the atom to test and its correct atomic
        number.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_position(case_data.atom, case_data.position)


def _test_get_position(atom, position):
    """
    Test :meth:`.PositionedAtom.get_position`

    Parameters
    ----------
    atom : :class:`.Atom`
        The atom to test.

    position : :class:`tuple` of :class:`float`
        The position (`x`, `y`, `z`) of the atom in cartesian
        coordinates.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert atom.get_position() == position
