import logging
from typing import Self

import stk

logger = logging.getLogger(__name__)


class PositionedAtom:
    """A container for stk.Atom and a coordinate.

    Parameters:
        atom:
            The atom.

        position:
            The position (`x`, `y`, `z`) of the atom in cartesian
            coordinates.

    """

    def __init__(
        self,
        atom: stk.Atom,
        position: tuple[float, ...],
    ) -> None:
        self._atom = atom
        self._position = position

    def get_atom(self) -> stk.Atom:
        return self._atom

    def get_atomic_number(self) -> int:
        return self._atom.get_atomic_number()

    def get_charge(self) -> int:
        return self._atom.get_charge()

    def get_id(self) -> int:
        return self._atom.get_id()

    def get_position(self) -> tuple[float, ...]:
        return self._position

    def _with_id(self, id: int) -> Self:  # noqa: A002
        """Modify the atom id."""
        self._atom = self._atom.with_id(id)
        return self

    def with_id(self, id: int) -> Self:  # noqa: A002
        """Get a clone but with a different id.

        Returns:
            A clone with a new id. Has the same type as the original
            atom.

        """
        return self.clone()._with_id(id)  # noqa: SLF001

    def clone(self) -> Self:
        """Return a clone.

        Returns:
            The clone. It has the same type as the original atom.

        """
        clone = self.__class__.__new__(self.__class__)
        clone._atom = self._atom  # noqa: SLF001
        clone._position = self._position  # noqa: SLF001
        return clone

    def __repr__(self) -> str:
        charge = (
            f", charge={self._atom.get_charge()}"
            if self._atom.get_charge() != 0
            else ""
        )
        return (
            f"{self._atom.__class__.__name__}({self._atom.get_id()}"
            f"{charge})"
        )

    def __str__(self) -> str:
        return repr(self)
