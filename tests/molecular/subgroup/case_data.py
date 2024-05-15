from dataclasses import dataclass

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    sub_group_data: dict[str, list[float]]
    name: str
