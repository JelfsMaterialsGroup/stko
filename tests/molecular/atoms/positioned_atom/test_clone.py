import stko
from tests.molecular.atoms.utilities import is_equivalent_atom


def test_clone(positioned_atom: stko.PositionedAtom) -> None:
    clone = positioned_atom.clone()
    is_equivalent_atom(positioned_atom, clone)
