"""
Torsion
=======

Class defining torsion angles.

"""

import logging

logger = logging.getLogger(__name__)


class Torsion:
    """
    Represents a torsion angle in a molecule.

    """

    def __init__(self, atom1, atom2, atom3, atom4):
        """
        Defines a torsion.

        Parameters
        ----------
        atom1 : :class:`stk.Atom`
            First atom in torsion.

        atom2 : :class:`stk.Atom`
            Second atom in torsion.

        atom3 : :class:`stk.Atom`
            Third atom in torsion.

        atom4 : :class:`stk.Atom`
            Fourth atom in torsion.

        """

        self._atom1 = atom1
        self._atom2 = atom2
        self._atom3 = atom3
        self._atom4 = atom4

    def get_atoms(self):
        return tuple((
            self._atom1,
            self._atom2,
            self._atom3,
            self._atom4,
        ))

    def get_atom_ids(self):
        return tuple(atom.get_id() for atom in self.get_atoms())

    def __iter__(self):
        return iter(self.get_atoms())

    def __str__(self):
        return (
            f'{self.__class__.__name__}({self._atom1}, '
            f'{self._atom2}, {self._atom3}, {self._atom4})'
        )

    def __repr__(self):
        return str(self)
