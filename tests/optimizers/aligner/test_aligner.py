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
    case_molecule.molecule.write('in1.mol')
    case_molecule.initial_molecule.write('in2.mol')
    optimizer = stko.Aligner(
        initial_molecule=case_molecule.initial_molecule,
        matching_pairs=(('C', 'C'), ('N', 'N')),
    )
    opt_res = optimizer.optimize(case_molecule.molecule)
    opt_res.write('out.mol')
    is_equivalent_molecule(opt_res, case_molecule.molecule)
    inequivalent_position_matrices(opt_res, case_molecule.molecule)

    test_rmsd = calculator.get_results(opt_res).get_rmsd()

    assert np.isclose(
        test_rmsd, case_molecule.rmsd, atol=1E-6
    )
    assert test_rmsd < test_rmsd_unopt
