import logging
from typing import Self

import stk

logger = logging.getLogger(__name__)


class Du:
    """Dummy of an stk.Atom.

    Parameters:
        id:
            ID of dummy atom.

    """

    def __init__(self, id: int) -> None:  # noqa: A002
        self._stk_atom = stk.Atom(
            id=id,
            atomic_number=1,
            charge=0,
        )
        self._atomic_number = 0
        self._charge = 0

    def get_id(self) -> int:
        return self._stk_atom.get_id()

    def get_atomic_number(self) -> int:
        return self._atomic_number

    def get_charge(self) -> int:
        return self._charge

    def clone(self) -> Self:
        return type(self)(self._stk_atom.get_id())

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._stk_atom.get_id()})"

    def __str__(self) -> str:
        return repr(self)
