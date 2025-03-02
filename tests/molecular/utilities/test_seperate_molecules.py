import stko

from .conftest import CaseData


def test_seperate_molecules(together_case_data: CaseData) -> None:
    distinct = stko.molecular_utilities.separate_molecule(
        together_case_data.molecule
    )
    assert together_case_data.molecule.get_num_atoms() == sum(
        [i[0].get_num_atoms() for i in distinct]
    )
    assert together_case_data.molecule.get_num_bonds() == sum(
        [i[0].get_num_bonds() for i in distinct]
    )
