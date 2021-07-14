import stko
from .utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)


def test_MMFF_opt(unoptimized_mol):
    optimizer = stko.MMFF()
    opt_res = optimizer.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)


def test_UFF_opt(unoptimized_mol):
    optimizer = stko.UFF()
    opt_res = optimizer.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)


def test_ETKDG_opt(unoptimized_mol):
    optimizer = stko.ETKDG()
    opt_res = optimizer.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)

    optimizer = stko.ETKDG(random_seed=2584)
    opt_res = optimizer.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)
