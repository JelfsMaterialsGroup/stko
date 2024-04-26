from dataclasses import dataclass

import stko


@dataclass(frozen=True, slots=True)
class CaseData:
    atom: stko.Du
    id: int
