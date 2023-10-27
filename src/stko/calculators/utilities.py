from collections import defaultdict

import numpy as np


def get_plane_normal(points):
    centroid = points.sum(axis=0) / len(points)
    return np.linalg.svd(points - centroid)[-1][2, :]


def get_atom_maps(mol):
    """
    Get atom maps from building blocks to constructude molecule.

    Returns a dictionary of dictionaries from atom id (in building
    block) to constructed molecule atom, indexed by building block id.

    Parameters
    ----------
    mol : :class:`.ConstructedMolecule`
        The :class:`.ConstructedMolecule` for which atom maps are
        desired.

    """
    atom_maps = defaultdict(dict)
    for atom_info in mol.get_atom_infos():
        bb_atom_id = atom_info.get_building_block_atom().get_id()
        atom_maps[atom_info.get_building_block_id()][
            bb_atom_id
        ] = atom_info.get_atom()
    return atom_maps


def is_inequivalent_atom(atom1, atom2):
    if atom1.__class__ is not atom2.__class__:
        return True
    if atom1.get_id() != atom2.get_id():
        return True
    if atom1.get_charge() != atom2.get_charge():
        return True
    if atom1.get_atomic_number() != atom2.get_atomic_number():
        return True
