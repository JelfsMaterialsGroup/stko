from scipy.spatial.distance import euclidean
import rdkit.Chem.AllChem as rdkit
import numpy as np
import itertools as it


def compute_angle(atom_a, atom_b, atom_c):
    # atom_a/b/c are numpy arrays for atomic positions.
    ba = atom_a - atom_b
    bc = atom_c - atom_b

    cosine_angle = np.dot(ba, bc) / (
        np.linalg.norm(ba) * np.linalg.norm(bc)
    )
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)


def compare_benzenes(initial_molecule, optimized_molecule):
    # Check position matrices.
    original_pos_mat = initial_molecule.get_position_matrix()
    new_pos_mat = optimized_molecule.get_position_matrix()
    assert not np.allclose(new_pos_mat, original_pos_mat)

    rdkit_mol = optimized_molecule.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)

    # Check C-C bond lengths.
    for atom_ids in rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts('[#6]~[#6]'),
    ):
        atoms = tuple(optimized_molecule.get_atoms(atom_ids))
        bond_length = euclidean(
            u=new_pos_mat[atoms[0].get_id()],
            v=new_pos_mat[atoms[1].get_id()],
        )

        assert abs(1.4 - bond_length) < 0.05

    # Check C-C-C bond angles
    for atom_ids in rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts('[#6]~[#6]~[#6]'),
    ):
        atoms = tuple(optimized_molecule.get_atoms(atom_ids))
        atom1_pos = new_pos_mat[atoms[0].get_id()]
        atom2_pos = new_pos_mat[atoms[1].get_id()]
        atom3_pos = new_pos_mat[atoms[2].get_id()]
        assert abs(
            compute_angle(atom1_pos, atom2_pos, atom3_pos) - 120
        ) < 5


def is_equivalent_atom(atom1, atom2):
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
    assert atom1.get_atomic_number() == atom2.get_atomic_number()
    assert atom1.__class__ is atom2.__class__


def is_equivalent_bond(bond1, bond2):
    assert bond1.__class__ is bond2.__class__
    assert bond1.get_order() == bond2.get_order()
    assert bond1.get_periodicity() == bond2.get_periodicity()
    is_equivalent_atom(bond1.get_atom1(), bond2.get_atom1())
    is_equivalent_atom(bond1.get_atom2(), bond2.get_atom2())


def is_equivalent_molecule(molecule1, molecule2):
    atoms = it.zip_longest(
        molecule1.get_atoms(),
        molecule2.get_atoms(),
    )
    for atom1, atom2 in atoms:
        is_equivalent_atom(atom1, atom2)

    bonds = it.zip_longest(
        molecule1.get_bonds(),
        molecule2.get_bonds(),
    )
    for bond1, bond2 in bonds:
        is_equivalent_bond(bond1, bond2)


def inequivalent_position_matrices(molecule1, molecule2):
    pos1 = molecule1.get_position_matrix()
    pos2 = molecule2.get_position_matrix()

    assert not np.allclose(pos1, pos2)
