"""
Fucntions for molecule conversion.

"""

import logging
import numpy as np
from ase import Atoms, Atom
from ase.io import write

from ...utilities import periodic_table

logger = logging.getLogger(__name__)


def stk_to_coord(mol, filename, cell=None):
    """
    Save stk.Molecule to a `.coord` file.

    """

    file_lines = []

    if cell is not None:
        file_lines.append('$periodic 3\n')
        # Write cell information.
        file_lines.append(
            '$cell angs\n'
            f'    {round(cell.a, 3)} {round(cell.b, 3)} '
            f'{round(cell.c, 3)} {round(cell.alpha, 3)} '
            f'{round(cell.beta, 3)} {round(cell.gamma, 3)}\n'
        )
    else:
        file_lines.append('')
        file_lines.append('')

    # Write coordinate information.
    file_lines.append('$coord angs\n')
    coords = mol.get_position_matrix()
    for atom in mol.get_atoms():
        coord = coords[atom.get_id()]
        element = atom.__class__.__name__

        file_lines.append(
            f'    {round(coord[0], 6)} {round(coord[1], 6)} '
            f'{round(coord[2], 6)} {element}\n'
        )

    file_lines.append('$end\n')

    with open(filename, 'w') as f:
        for line in file_lines:
            f.write(line)


def stk_to_ase(mol, cell=None):
    """
    Return ase.Atoms from stk.Molecule.

    """

    ase_atoms = Atoms()

    coords = mol.get_position_matrix()
    for atom in mol.get_atoms():
        coord = coords[atom.get_id()]
        element = atom.__class__.__name__
        ase_atoms.append(Atom(symbol=element, position=coord))

    if cell is not None:
        ase_atoms.set_cell([
            cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma
        ])
        ase_atoms.set_pbc(True)

    return ase_atoms


def write_cif(mol, cell, filename):
    """
    Write CIF from stk.Molecule and .Cell in P1 symmetry.

    This function uses the ASE interface to write a CIF.

    """

    # Convert stk Molecule object into ase Atoms object.
    ase_atoms = stk_to_ase(mol, cell)

    write(images=ase_atoms, filename=filename)


def with_structure_from_periodic_turbomole(mol, path):
    """
    Return a clone, with its structure taken from a Turbomole file.

    Note that coordinates in ``.coord`` files are given in Bohr.

    Parameters
    ----------
    path : :class:`str`
        The full path of the ``.coord`` file from which the
        structure should be updated.

    Returns
    -------
    :class:`.Molecule`
        A clone with atomic positions found in `path`.

    Raises
    ------
    :class:`RuntimeError`
        If the number of atoms in the file does not match the
        number of atoms in the molecule or if atom elements in the
        file do not agree with the atom elements in the molecule.

    """

    bohr_to_ang = 0.5291772105638411

    content = []
    with open(path, 'r') as f:
        for line in f.readlines()[1:-6]:
            content.append(line)

    # Check the atom count is correct.
    num_atoms = mol.get_num_atoms()
    if len(content) != num_atoms:
        raise RuntimeError(
            'The number of atoms in the coord file, '
            f'{len(content)}, does not match the number of atoms '
            f'in the molecule, {num_atoms}.'
        )

    # Save all the coords in the file.
    new_coords = []
    elements = [i.__class__.__name__ for i in mol.get_atoms()]
    for i, line in enumerate(content):
        *coords, element = line.split()
        if element.isnumeric():
            element = periodic_table[int(element)]

        if element != elements[i]:
            raise RuntimeError(
                f'Atom {i} element does not match file.'
            )

        new_coords.append([float(i)*bohr_to_ang for i in coords])

    # Check that the correct number of atom
    # lines was present in the file.
    if i+1 != num_atoms:
        raise RuntimeError(
            'The number of atoms lines in the coord file, '
            f'{i+1}, does not match the number of atoms '
            f'in the molecule, {num_atoms}.'
        )

    # Update the structure.
    return mol.with_position_matrix(np.asarray(new_coords))
