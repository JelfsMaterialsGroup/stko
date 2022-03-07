from turtle import width
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
    optimizer = stko.Aligner(
        initial_molecule=case_molecule.initial_molecule,
        matching_pairs=(('C', 'C'), ('N', 'N')),
    )
    opt_res = optimizer.optimize(case_molecule.molecule)
    is_equivalent_molecule(opt_res, case_molecule.molecule)
    inequivalent_position_matrices(opt_res, case_molecule.molecule)

    test_rmsd = calculator.get_results(opt_res).get_rmsd()

    assert np.isclose(
        test_rmsd, case_molecule.rmsd, atol=1E-6
    )
    assert test_rmsd < test_rmsd_unopt


def test_alignment_potential(case_potential):
    aligner = stko.Aligner(
        initial_molecule=case_potential.initial_molecule,
        matching_pairs=case_potential.pairs,
    )
    potential = stko.AlignmentPotential(
        matching_pairs=case_potential.pairs,
        width=2,
    )
    supramolecule = aligner._get_supramolecule(case_potential.molecule)
    assert case_potential.potential == potential.compute_potential(
        supramolecule,
    )
