from dataclasses import dataclass

import stko


@dataclass(frozen=True, slots=True)
class CaseData:
    atom: stko.PositionedAtom
    id: int
    charge: int
    atomic_number: int
    position: tuple[float, float, float]
