import numpy as np
import stk

import stko
from tests.optimizers.rdkit.conftest import CaseData
from tests.optimizers.utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)


def test_mmff_opt(case_mmff_molecule: CaseData) -> None:
    energy = stko.MMFFEnergy().get_energy(case_mmff_molecule.molecule)
    optimizer = stko.MMFF()
    opt_res = optimizer.optimize(case_mmff_molecule.molecule)
    is_equivalent_molecule(opt_res, case_mmff_molecule.molecule)
    inequivalent_position_matrices(opt_res, case_mmff_molecule.molecule)
    opt_energy = stko.MMFFEnergy().get_energy(opt_res)
    assert np.isclose(energy, case_mmff_molecule.unoptimised_energy, atol=1e-3)
    assert opt_energy < case_mmff_molecule.unoptimised_energy


def test_uff_opt(case_uff_molecule: CaseData) -> None:
    energy = stko.UFFEnergy().get_energy(case_uff_molecule.molecule)
    optimizer = stko.UFF()
    opt_res = optimizer.optimize(case_uff_molecule.molecule)
    is_equivalent_molecule(opt_res, case_uff_molecule.molecule)
    inequivalent_position_matrices(opt_res, case_uff_molecule.molecule)
    opt_energy = stko.UFFEnergy().get_energy(opt_res)
    assert np.isclose(energy, case_uff_molecule.unoptimised_energy, atol=1e-3)
    assert opt_energy < case_uff_molecule.unoptimised_energy


def test_etkdg_opt(case_etkdg_molecule: stk.BuildingBlock) -> None:
    optimizer = stko.ETKDG()
    opt_res = optimizer.optimize(case_etkdg_molecule)
    is_equivalent_molecule(opt_res, case_etkdg_molecule)
    inequivalent_position_matrices(opt_res, case_etkdg_molecule)

    optimizer = stko.ETKDG(random_seed=2584)
    opt_res = optimizer.optimize(case_etkdg_molecule)
    is_equivalent_molecule(opt_res, case_etkdg_molecule)
    inequivalent_position_matrices(opt_res, case_etkdg_molecule)
