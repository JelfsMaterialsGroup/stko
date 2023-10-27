import numpy as np
import pytest
import stko


def test_rmsd(case_data):
    calculator = stko.RmsdCalculator(case_data.mol1)
    results = calculator.get_results(case_data.mol2)
    test_rmsd = results.get_rmsd()
    assert np.isclose(test_rmsd, case_data.rmsd, atol=1e-4)


def test_rmsd_ignore_hydrogens(ignore_h_case_data):
    calculator = stko.RmsdCalculator(ignore_h_case_data.mol1, True)
    results = calculator.get_results(ignore_h_case_data.mol2)
    test_rmsd = results.get_rmsd()
    assert np.isclose(test_rmsd, ignore_h_case_data.rmsd, atol=1e-4)


def test_rmsd_aligned(aligned_case_data):
    calculator = stko.RmsdMappedCalculator(aligned_case_data.mol1)
    results = calculator.get_results(aligned_case_data.mol2)
    test_rmsd = results.get_rmsd()
    assert np.isclose(test_rmsd, aligned_case_data.rmsd, atol=1e-4)


def test_rmsd_different_molecule(different_case_data):
    calculator = stko.RmsdCalculator(different_case_data.mol1)
    with pytest.raises(stko.DifferentMoleculeError):
        calculator.get_results(different_case_data.mol2)


def test_rmsd_different_atom(ordering_case_data):
    calculator = stko.RmsdCalculator(ordering_case_data.mol1)
    with pytest.raises(stko.DifferentAtomError):
        calculator.get_results(ordering_case_data.mol2)
