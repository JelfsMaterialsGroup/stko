import stko
import pytest

from ..utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)
import numpy as np

try:
    from openbabel import openbabel
except ImportError:
    openbabel = None


def test_open_babel(case_molecule):

    if openbabel is None:
        with pytest.raises(stko.WrapperNotInstalledException) as e:
            energy = (
                stko.OpenBabelEnergy('uff').get_energy(
                    mol=case_molecule.molecule,
                )
            )
    else:
        energy = (
            stko.OpenBabelEnergy('uff').get_energy(
                mol=case_molecule.molecule,
            )
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
