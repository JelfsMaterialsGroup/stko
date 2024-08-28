import stk

import stko
from tests.molecular.decompose.case_data import CaseData


def test_decomposemoc(case_data: CaseData) -> None:
    ligands = stko.molecule_analysis.DecomposeMOC().decompose(
        molecule=case_data.cage,
        metal_atom_nos=case_data.metal_atom_nos,
    )

    assert case_data.num_ligands == len(ligands)

    for ligand in ligands:
        smiles = stk.Smiles().get_key(ligand)
        assert smiles in case_data.bb_smiles
