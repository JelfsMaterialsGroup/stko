from ..utilities import is_equivalent_atom


def test_clone(dummy_atom):
    """
    Test :meth:`.Du.clone`.

    Parameters
    ----------
    atom : :class:`.PositionedAtom`
        The atom to be cloned.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = dummy_atom.clone()
    is_equivalent_atom(dummy_atom, clone)
