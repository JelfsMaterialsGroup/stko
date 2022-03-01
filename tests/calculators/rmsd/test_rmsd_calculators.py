import stko
import pytest
import numpy as np


def test_rmsd(case_data):
    calculator = stko.RmsdCalculator(case_data.mol1)
    results = calculator.get_results(case_data.mol2)
    test_rmsd = results.get_rmsd()
    assert np.isclose(test_rmsd, case_data.rmsd, atol=1E-4)


def test_rmsd_ignore_hydrogens(ignore_h_case_data):
    calculator = stko.RmsdCalculator(ignore_h_case_data.mol1, True)
    results = calculator.get_results(ignore_h_case_data.mol2)
    test_rmsd = results.get_rmsd()
    assert test_rmsd == ignore_h_case_data.rmsd


def test_rmsd_different_molecule(different_case_data):
    calculator = stko.RmsdCalculator(different_case_data.mol1)
    with pytest.raises(stko.DifferentMoleculeException):
        calculator.get_results(different_case_data.mol2)


def test_rmsd_different_atom(ordering_case_data):
    calculator = stko.RmsdCalculator(ordering_case_data.mol1)
    with pytest.raises(stko.DifferentAtomException):
        calculator.get_results(ordering_case_data.mol2)
