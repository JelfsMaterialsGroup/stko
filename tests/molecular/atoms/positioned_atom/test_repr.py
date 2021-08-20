import stk
import stko

from ..utilities import is_equivalent_atom


def test_repr(positioned_atom):
    """
    Test :meth:`.PositionedAtom.__repr__`.

    Parameters
    ----------
    atom : :class:`.PositionedAtom`
        The atom, whose representation should be tested.

    Returns
    -------
    None : :class:`NoneType`

    """

    other = stko.PositionedAtom(
        atom=eval(repr(positioned_atom), dict(stk.__dict__)),
        position=positioned_atom.get_position(),
    )
    is_equivalent_atom(other, positioned_atom)
