from scipy.spatial.distance import euclidean
import numpy as np


def compute_angle(atom_a, atom_b, atom_c):
    # atom_a/b/c are numpy arrays for atomic positions
    ba = atom_a - atom_b
    bc = atom_c - atom_b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)


def compare_molecules(initial_molecule, optimized_molecule):
    # Check position matrices.
    original_pos_mat = initial_molecule.get_position_matrix()
    new_pos_mat = optimized_molecule.get_position_matrix()
    assert not np.allclose(new_pos_mat, original_pos_mat)

    # Check C-C bond lengths.
    for bond in initial_molecule.get_bonds():
        atom1 = bond.get_atom1()
        atom2 = bond.get_atom2()
        if (
            atom1.get_atomic_number() == 6
            and atom2.get_atomic_number() == 6
        ):
            bond_length = euclidean(
                u=new_pos_mat[atom1.get_id()],
                v=new_pos_mat[atom2.get_id()],
            )

            assert abs(1.4 - bond_length) < 0.05
    
    # Check C-C-C bond angles
    angles = {x: [] for x in initial_molecule.get_atoms()}

    # Get angle triplets
    for bond in initial_molecule.get_bonds():
        atom1 = bond.get_atom1()
        atom2 = bond.get_atom2()
        angles[atom1].append(atom2)
        angles[atom2].append(atom1)
    
    for angle in angles.items():
        atom_a = new_pos_mat[angle[1][0].get_id()]
        atom_b = new_pos_mat[angle[0].get_id()]
        atom_c = new_pos_mat[angle[1][1].get_id()]
        assert abs(compute_angle(atom_a, atom_b, atom_c) - 120) < 5
