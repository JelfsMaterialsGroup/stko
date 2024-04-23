from typing import TypeVar

import stk

MoleculeT = TypeVar("MoleculeT", bound=stk.Molecule)
ConstructedMoleculeT = TypeVar(
    "ConstructedMoleculeT", bound=stk.ConstructedMolecule
)
