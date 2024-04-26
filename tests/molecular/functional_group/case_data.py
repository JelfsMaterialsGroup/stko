from dataclasses import dataclass

import stko


@dataclass(frozen=True, slots=True)
class CaseData:
    smiles: str
    num_fgs: int
    factory: stko.functional_groups.ThreeSiteFactory
    name: str
