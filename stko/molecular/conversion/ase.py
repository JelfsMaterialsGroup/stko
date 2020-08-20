"""
Functions for molecule conversion between stk and ASE.

"""

import logging
from ase import Atoms, Atom

logger = logging.getLogger(__name__)


class ASEConverter:
    """
    Converters for the atomic simulation environment (ASE).

    """

    def to_stk(self, atoms, periodic=False):
        """
        Convert ase.Atoms() to stk.Molecule.

        """

        raise NotImplementedError()

    def from_stk(self, molecule, unit_cell=None):
        """
        Convert molecule to ase.Atoms().

        """

        ase_atoms = Atoms()

        coords = molecule.get_position_matrix()
        for atom in molecule.get_atoms():
            coord = coords[atom.get_id()]
            element = atom.__class__.__name__
            ase_atoms.append(Atom(symbol=element, position=coord))

        if unit_cell is not None:
            ase_atoms.set_cell([
                unit_cell.get_a(),
                unit_cell.get_b(),
                unit_cell.get_c(),
                unit_cell.get_alpha(),
                unit_cell.get_beta(),
                unit_cell.get_gamma(),
            ])
            ase_atoms.set_pbc(True)

        return ase_atoms
