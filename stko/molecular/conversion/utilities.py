"""
Utilities for molecule conversion.

"""

import logging
import numpy as np
from ase.io import write
from .ase import ASEConverter

from ...utilities import periodic_table

logger = logging.getLogger(__name__)


def write_cif(mol, unit_cell, filename):
    """
    Write CIF from stk.Molecule and .Cell in P1 symmetry.

    This function uses the ASE interface to write a CIF.

    """

    # Convert stk Molecule object into ase Atoms object.
    ase_atoms = ASEConverter().from_stk(
        molecule=mol,
        unit_cell=unit_cell,
    )

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

    :class:`RuntimeError`
        If the the Turbomole file has coordinates defined as fractional
        based on a unit cell.

    """

    bohr_to_ang = 0.5291772105638411
    num_atoms = mol.get_num_atoms()
    with open(path, 'r') as f:
        content = f.readlines()

    periodicity = False
    lattice_vectors = None
    lattice_units = None
    cell_parameters = None
    cell_units = None
    for line_number, line in enumerate(content):
        if '$periodic' in line:
            periodicity = int(line.rstrip().split()[1])
        if '$cell' in line:
            if 'angs' in line:
                cell_units = 'angstrom'
            else:
                raise ValueError('cell not in Angstroms.')
            cell_parameters = [
                float(j)
                for j in content[line_number+1].rstrip().split()
            ]
        if '$lattice' in line:
            if 'angs' in line:
                lattice_units = 'angstrom'
            else:
                raise ValueError('lattice not in Angstroms.')
            lattice_vectors = (
                [
                    float(j)
                    for j in content[line_number+1].rstrip().split()
                ],
                [
                    float(j)
                    for j in content[line_number+2].rstrip().split()
                ],
                [
                    float(j)
                    for j in content[line_number+3].rstrip().split()
                ],
            )
        if '$coord' in line:
            if 'angs' in line:
                coord_units = 'angstrom'
            elif 'frac' in line:
                coord_units = 'fractional'
                raise RuntimeError(
                    'Fractional coordinates are not handled currently.'
                )
            elif 'bohr' in line:
                coord_units = 'bohr'
            else:
                coord_units = 'bohr'
            coord_section = (
                content[line_number+1:line_number+1+num_atoms]
            )

    # Check that cell is only defined once.
    chk2 = lattice_vectors is not None and cell_parameters is not None
    if periodicity and chk2:
        raise RuntimeError('The cell is defined twice in the file.')

    # Save all the coords in the file.
    new_coords = []
    elements = [i.__class__.__name__ for i in mol.get_atoms()]
    for _id, line in enumerate(coord_section):
        *coords, element = line.split()
        if element.isnumeric():
            element = periodic_table[int(element)]

        if element != elements[i]:
            raise RuntimeError(
                f'Atom {_id} element does not match file.'
            )
        if coord_units == 'bohr':
            new_coords.append([float(i)*bohr_to_ang for i in coords])
        elif coord_units == 'fractional':
            raise RuntimeError(
                'Fractional coordinates are not handled currently.'
            )
        else:
            new_coords.append([float(i) for i in coords])

    # Check that the correct number of atom
    # lines was present in the file.
    if _id+1 != num_atoms:
        raise RuntimeError(
            'The number of atoms lines in the coord file, '
            f'{_id+1}, does not match the number of atoms '
            f'in the molecule, {num_atoms}.'
        )

    # Update the structure.
    return mol.with_position_matrix(np.asarray(new_coords))
