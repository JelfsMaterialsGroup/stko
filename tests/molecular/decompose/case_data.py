from dataclasses import dataclass

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    cage: stk.Molecule
    num_ligands: int
    metal_atom_nos: tuple[int]
    bb_smiles: tuple[str]
    name: str
