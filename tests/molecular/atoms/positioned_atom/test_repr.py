import stk

import stko
from tests.molecular.atoms.utilities import is_equivalent_atom


def test_repr(positioned_atom: stko.PositionedAtom) -> None:
    other = stko.PositionedAtom(
        atom=eval(repr(positioned_atom), dict(stk.__dict__)),  # noqa: S307
        position=positioned_atom.get_position(),
    )
    is_equivalent_atom(other, positioned_atom)
