import numpy as np

import stko
from tests.optimizers.aligner.conftest import CaseData, CasePotential
from tests.optimizers.utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)


def test_aligner(case_molecule: CaseData) -> None:
    calculator = stko.RmsdCalculator(case_molecule.initial_molecule)
    test_rmsd_unopt = calculator.get_results(case_molecule.molecule).get_rmsd()
    optimizer = stko.Aligner(
        initial_molecule=case_molecule.initial_molecule,
        matching_pairs=(("C", "C"), ("N", "N")),
    )
    opt_res = optimizer.optimize(case_molecule.molecule)
    is_equivalent_molecule(opt_res, case_molecule.molecule)
    inequivalent_position_matrices(opt_res, case_molecule.molecule)

    test_rmsd = calculator.get_results(opt_res).get_rmsd()

    assert np.isclose(test_rmsd, case_molecule.rmsd, atol=1e-6)
    assert test_rmsd < test_rmsd_unopt


def test_alignment_potential(case_potential: CasePotential) -> None:
    aligner = stko.Aligner(
        initial_molecule=case_potential.initial_molecule,
        matching_pairs=case_potential.pairs,
    )
    potential = stko.AlignmentPotential(
        matching_pairs=case_potential.pairs,
        width=2,
    )
    supramolecule = aligner._get_supramolecule(case_potential.molecule)  # noqa: SLF001
    assert np.isclose(
        case_potential.potential,
        potential.compute_potential(
            supramolecule,
        ),
        atol=1e-6,
    )
