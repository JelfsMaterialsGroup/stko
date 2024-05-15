import stko

from .case_data import CaseData


def test_get_position(case_data: CaseData) -> None:
    _test_get_position(case_data.atom, case_data.position)


def _test_get_position(
    atom: stko.PositionedAtom, position: tuple[float, float, float]
) -> None:
    assert atom.get_position() == position
