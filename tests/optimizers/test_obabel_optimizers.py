import stko
from .utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)


def test_open_babel(unoptimized_mol):
    optimizer = stko.OpenBabel('uff')
    opt_res = optimizer.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)
