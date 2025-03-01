"""Module of molecular utilities."""

from collections import abc

import numpy as np
import stk
from rdkit.Chem import AllChem as rdkit  # noqa: N813


def merge_stk_molecules(
    molecules: abc.Sequence[stk.Molecule],
) -> stk.BuildingBlock:
    """Merge any list of separate molecules to be in one class."""
    atoms: list[stk.Atom] = []
    bonds: list[stk.Bond] = []
    pos_mat = []
    for molecule in molecules:
        atom_ids_map = {}
        for atom in molecule.get_atoms():
            new_id = len(atoms)
            atom_ids_map[atom.get_id()] = new_id
            atoms.append(
                stk.Atom(
                    id=atom_ids_map[atom.get_id()],
                    atomic_number=atom.get_atomic_number(),
                    charge=atom.get_charge(),
                )
            )

        bonds.extend(
            i.with_ids(id_map=atom_ids_map) for i in molecule.get_bonds()
        )
        pos_mat.extend(list(molecule.get_position_matrix()))

    return stk.BuildingBlock.init(
        atoms=atoms,
        bonds=bonds,
        position_matrix=np.array(pos_mat),
    )


def update_stk_from_rdkit_conformer(
    stk_mol: stk.Molecule,
    rdk_mol: rdkit.Mol,
    conf_id: int,
) -> stk.Molecule:
    """Update the structure to match `conf_id` of `mol`.

    Parameters:
        struct:
            The molecule whoce coordinates are to be updated.

        mol:
            The :mod:`rdkit` molecule to use for the structure update.

        conf_id:
            The conformer ID of the `mol` to update from.

    Returns:
        The molecule.

    TODO: Add to stko.

    """
    pos_mat = rdk_mol.GetConformer(id=conf_id).GetPositions()
    return stk_mol.with_position_matrix(pos_mat)
