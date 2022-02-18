import stko
from ..utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)
import numpy as np


def test_open_babel(case_molecule):
    energy = (
        stko.OpenBabelEnergy('uff').get_energy(case_molecule.molecule)
    )
    optimizer = stko.OpenBabel('uff')
    opt_res = optimizer.optimize(case_molecule.molecule)
    is_equivalent_molecule(opt_res, case_molecule.molecule)
    inequivalent_position_matrices(opt_res, case_molecule.molecule)
    opt_energy = (
        stko.OpenBabelEnergy('uff').get_energy(opt_res)
    )
    assert np.isclose(
        energy, case_molecule.unoptimised_energy, atol=1E-3
    )
    assert opt_energy < case_molecule.unoptimised_energy
