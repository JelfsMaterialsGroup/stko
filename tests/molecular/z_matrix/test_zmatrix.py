import stk

import stko
from tests.molecular.z_matrix.case_data import CaseData


def test_zmatrix(case_data: CaseData) -> None:
    _test_zmatrix(case_data.molecule, case_data.zmatrix)


def _test_zmatrix(molecule: stk.Molecule, zmatrix: str) -> None:
    result = stko.ZMatrix().get_zmatrix(molecule)
    assert result == zmatrix
