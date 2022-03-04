import stko
from ..utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)
import numpy as np


def test_aligner(case_molecule):
    calculator = stko.RmsdCalculator(case_molecule.initial_molecule)
    test_rmsd_unopt = calculator.get_results(
        case_molecule.molecule
    ).get_rmsd()

    optimizer = stko.Aligner(case_molecule.initial_molecule)
    opt_res = optimizer.optimize(case_molecule.molecule)
    is_equivalent_molecule(opt_res, case_molecule.molecule)
    inequivalent_position_matrices(opt_res, case_molecule.molecule)

    test_rmsd = calculator.get_results(opt_res).get_rmsd()

    assert np.isclose(
        test_rmsd, case_molecule.rmsd, atol=1E-6
    )
    assert test_rmsd < test_rmsd_unopt
