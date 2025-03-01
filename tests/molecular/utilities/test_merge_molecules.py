
import stko
from .conftest import CaseData


def test_merge_molecules(case_data: CaseData) -> None:
    merged = stko.molecular_utilities.merge_stk_molecules(case_data.molecules)
    assert merged.get_num_atoms() == sum(
        [i.get_num_atoms() for i in case_data.molecules]
    )
    assert merged.get_num_bonds() == sum(
        [i.get_num_bonds() for i in case_data.molecules]
    )
