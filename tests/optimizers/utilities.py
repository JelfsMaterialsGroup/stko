from scipy.spatial.distance import euclidean
import numpy as np


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
