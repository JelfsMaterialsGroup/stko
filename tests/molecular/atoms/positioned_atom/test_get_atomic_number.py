import stko

from .case_data import CaseData


def test_get_atomic_number(case_data: CaseData) -> None:
    _test_get_atomic_number(case_data.atom, case_data.atomic_number)


def _test_get_atomic_number(
    atom: stko.PositionedAtom, atomic_number: int
) -> None:
    assert atom.get_atomic_number() == atomic_number
