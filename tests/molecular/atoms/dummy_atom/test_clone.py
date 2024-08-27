import stko
from tests.molecular.atoms.utilities import is_equivalent_atom


def test_clone(dummy_atom: stko.Du) -> None:
    clone = dummy_atom.clone()
    is_equivalent_atom(dummy_atom, clone)
