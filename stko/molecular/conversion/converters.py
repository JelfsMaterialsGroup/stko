"""
Fucntions for molecule conversion.

"""

import logging
from ase import Atoms, Atom
from ase.io import write

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
