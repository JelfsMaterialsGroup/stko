import stko
from .utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)


def test_open_babel(
    unoptimized_mol, unoptimized_obabel_uff, optimized_obabel_uff,
):
    energy = (
        stko.OpenBabelEnergy('uff').get_energy(unoptimized_mol)
    )
    optimizer = stko.OpenBabel('uff')
    opt_res = optimizer.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)
    opt_energy = (
        stko.OpenBabelEnergy('uff').get_energy(opt_res)
    )
    assert energy == unoptimized_obabel_uff
    assert opt_energy == optimized_obabel_uff
