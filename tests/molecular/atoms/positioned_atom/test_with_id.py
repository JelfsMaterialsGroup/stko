import stko
from tests.molecular.atoms.utilities import is_equivalent_atom


def test_with_id(positioned_atom: stko.PositionedAtom, id: int) -> None:  # noqa: A002
    # Save a clone to ensure that "atom" is not changed by the test.
    before = positioned_atom.clone()
    _test_with_id(positioned_atom, id)
    is_equivalent_atom(before, positioned_atom)


def _test_with_id(positioned_atom: stko.PositionedAtom, id: int) -> None:  # noqa: A002
    new_atom = positioned_atom.with_id(id)
    assert new_atom is not positioned_atom
    assert new_atom.__class__ is positioned_atom.__class__
    assert new_atom.get_id() == id
    assert new_atom.get_charge() == positioned_atom.get_charge()
    atm_num = positioned_atom.get_atomic_number()
    assert new_atom.get_atomic_number() == atm_num
    assert new_atom.get_position() == positioned_atom.get_position()
