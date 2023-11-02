import logging
import typing

import stk

logger = logging.getLogger(__name__)


class Torsion:
    """
    Represents a torsion angle in a molecule.

    """

    def __init__(
        self,
        atom1: stk.Atom,
        atom2: stk.Atom,
        atom3: stk.Atom,
        atom4: stk.Atom,
    ) -> None:
        """
        Defines a torsion.

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

        self._atom1 = atom1
        self._atom2 = atom2
        self._atom3 = atom3
        self._atom4 = atom4

    def get_atoms(self) -> tuple[stk.Atom, ...]:
        return tuple(
            (
                self._atom1,
                self._atom2,
                self._atom3,
                self._atom4,
            )
        )

    def get_atom_ids(self) -> tuple[int, ...]:
        return tuple(atom.get_id() for atom in self.get_atoms())

    def __iter__(self) -> typing.Iterable[stk.Atom]:
        return iter(self.get_atoms())

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}({self._atom1}, "
            f"{self._atom2}, {self._atom3}, {self._atom4})"
        )

    def __repr__(self) -> str:
        return str(self)
