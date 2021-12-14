from ..utilities import is_equivalent_atom


def test_clone(positioned_atom):
    """
    Test :meth:`.PositionedAtom.clone`.

    Parameters
    ----------
    atom : :class:`.PositionedAtom`
        The atom to be cloned.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = positioned_atom.clone()
    is_equivalent_atom(positioned_atom, clone)
