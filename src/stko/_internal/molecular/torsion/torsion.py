import logging
from collections import abc

import stk

logger = logging.getLogger(__name__)


class Torsion:
    """Represents a torsion angle in a molecule.

    Parameters:
        atom1:
            First atom in torsion.

        atom2:
            Second atom in torsion.

        atom3:
            Third atom in torsion.

        atom4:
            Fourth atom in torsion.

    """

    def __init__(
        self,
        atom1: stk.Atom | None,
        atom2: stk.Atom | None,
        atom3: stk.Atom | None,
        atom4: stk.Atom | None,
    ) -> None:
        self._atom1 = atom1
        self._atom2 = atom2
        self._atom3 = atom3
        self._atom4 = atom4

    def get_atoms(self) -> tuple[stk.Atom | None, ...]:
        return (
            self._atom1,
            self._atom2,
            self._atom3,
            self._atom4,
        )

    def get_atom_ids(self) -> tuple[int, ...]:
        return tuple(
            atom.get_id()  # type: ignore[union-attr]
            for atom in self.get_atoms()
        )

    def __iter__(self) -> abc.Iterator[stk.Atom | None]:
        return iter(self.get_atoms())

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}({self._atom1}, "
            f"{self._atom2}, {self._atom3}, {self._atom4})"
        )

    def __repr__(self) -> str:
        return str(self)
