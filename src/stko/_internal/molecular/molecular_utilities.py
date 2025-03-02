"""Module of molecular utilities."""

from collections import abc

import numpy as np
import stk
from rdkit.Chem import AllChem as rdkit  # noqa: N813

from stko._internal.molecular.networkx.network import Network


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


def separate_molecule(
    molecule: stk.Molecule,
) -> abc.Sequence[tuple[stk.BuildingBlock, list[int]]]:
    """Given a molecule, it returns distinct disconnected molecules."""
    network = Network.init_from_molecule(molecule)
    connected = network.get_connected_components()
    molecules = []
    for cg in connected:
        # Get atoms from nodes.
        atoms = list(cg)
        atom_ids = tuple(i.get_id() for i in atoms)

        # Sort both by atom id.
        atom_ids, atoms = zip(  # type: ignore[assignment]
            *sorted(zip(atom_ids, atoms, strict=True)), strict=True
        )

        atom_ids_map = {atom_ids[i]: i for i in range(len(atom_ids))}
        new_mol = stk.BuildingBlock.init(
            atoms=(
                stk.Atom(
                    id=atom_ids_map[i.get_id()],
                    atomic_number=i.get_atomic_number(),
                    charge=i.get_charge(),
                )
                for i in atoms
            ),
            bonds=(
                i.with_ids(id_map=atom_ids_map)
                for i in molecule.get_bonds()
                if i.get_atom1().get_id() in atom_ids
                and i.get_atom2().get_id() in atom_ids
            ),
            position_matrix=np.array(
                tuple(
                    i for i in molecule.get_atomic_positions(atom_ids=atom_ids)
                )
            ),
        )
        molecules.append((new_mol, list(atom_ids)))
    return molecules


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
