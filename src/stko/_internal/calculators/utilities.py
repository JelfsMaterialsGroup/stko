from collections import defaultdict

import numpy as np
import stk


def get_plane_normal(points: np.ndarray) -> np.ndarray:
    centroid = points.sum(axis=0) / len(points)
    return np.linalg.svd(points - centroid)[-1][2, :]


def get_atom_maps(mol: stk.ConstructedMolecule) -> dict:
    """Get atom maps from building blocks to constructude molecule.

    Returns a dictionary of dictionaries from atom id (in building
    block) to constructed molecule atom, indexed by building block id.

    Parameters:
        mol:
            The :class:`stk.ConstructedMolecule` for which atom maps are
            desired.

    """
    atom_maps = defaultdict(dict)  # type: ignore[var-annotated]
    for atom_info in mol.get_atom_infos():
        if atom_info is not None:
            bb_atom = atom_info.get_building_block_atom()
            bb_atom_id = bb_atom.get_id()  # type: ignore[union-attr]
            atom_maps[atom_info.get_building_block_id()][bb_atom_id] = (
                atom_info.get_atom()
            )
    return atom_maps


def is_inequivalent_atom(atom1: stk.Atom, atom2: stk.Atom) -> bool | None:
    if atom1.__class__ is not atom2.__class__:
        return True
    if atom1.get_id() != atom2.get_id():
        return True
    if atom1.get_charge() != atom2.get_charge():
        return True
    if atom1.get_atomic_number() != atom2.get_atomic_number():
        return True
    return None
