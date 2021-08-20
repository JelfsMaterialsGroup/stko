"""
Dummy Atom
==========

#. :class:`.Du`

A class representing a dummy stk.Atom.

"""

import logging
import stk

logger = logging.getLogger(__name__)


class Du:
    """
    Dummy of an stk.Atom.

    """

    def __init__(self, id):
        """
        Initialize a :class:`.Du`

        Parameters
        ----------

        """

        self._stk_atom = stk.Atom(
            id=id,
            atomic_number=1,
            charge=0,
        )
        self._atomic_number = 0
        self._charge = 0

    def get_id(self):
        """
        Get the ID of the atom.

        Returns
        ------
        :class:`int`
            The atom ID.

        """

        return self._stk_atom.get_id()

    def get_atomic_number(self):
        """
        Get the atomic number of the atom.

        Returns
        ------
        :class:`int`
            The atomic number.

        """

        return self._atomic_number

    def get_charge(self):
        """
        Get the charge of the atom.

        Returns
        ------
        :class:`int`
            The atom charge.

        """

        return self._charge

    def clone(self):
        """
        Returns a clone.

        Returns
        ------
        :class:`.Atom`
            A clone.

        """

        return type(self)(self._stk_atom.get_id())

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self._stk_atom.get_id()})'

    def __str__(self) -> str:
        return repr(self)
