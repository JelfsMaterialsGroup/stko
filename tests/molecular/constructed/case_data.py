from dataclasses import dataclass

import numpy as np
import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    constructed_molecule: stk.ConstructedMolecule
    centroids: dict[int, np.ndarray]
    atom_ids: dict[int, list[int]]
    name: str
