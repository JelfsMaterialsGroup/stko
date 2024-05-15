from dataclasses import dataclass

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    zmatrix: str
