"""
Positioned Atom
===============

#. :class:`.PositionedAtom`

A class representing a positioned stk.Atom.

"""

import logging

logger = logging.getLogger(__name__)


class PositionedAtom:
    """
    A container for stk.Atom and a coordinate.

    """

    def __init__(self, atom, position):
        """
        Initialize a :class:`PositionedAtom`.

        Parameters
        ----------
        atom : :class:`stk.Atom`
            The atom.

        position : :class:`tuple` of :class:`float`
            The position (`x`, `y`, `z`) of the atom in cartesian
            coordinates.

        """

        self._atom = atom
        self._position = position

    def get_atom(self):
        return self._atom

    def get_atomic_number(self):
        return self._atom.get_atomic_number()

    def get_charge(self):
        return self._atom.get_charge()

    def get_id(self):
        return self._atom.get_id()

    def get_position(self):
        return self._position

    def _with_id(self, id):
        """
        Modify the atom id.

        """

        self._atom = self._atom.with_id(id)
        return self

    def with_id(self, id):
        """
        Get a clone but with a different id.

        Returns
        -------
        :class:`.Atom`
            A clone with a new id. Has the same type as the original
            atom.

        """

        return self.clone()._with_id(id)

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.Atom`
            The clone. It has the same type as the original atom.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._atom = self._atom
        clone._position = self._position
        return clone

    def __repr__(self):
        charge = (
            f', charge={self._atom.get_charge()}'
            if self._atom.get_charge() != 0 else ''
        )
        return (
            f'{self._atom.__class__.__name__}({self._atom.get_id()}'
            f'{charge})'
        )

    def __str__(self):
        return repr(self)
