from scipy.spatial.distance import euclidean
import rdkit.Chem.AllChem as rdkit
import numpy as np


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
