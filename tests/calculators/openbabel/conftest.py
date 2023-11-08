import pytest
import stk


def a_molecule():
    return stk.BuildingBlock(smiles="CCCCCC")


@pytest.fixture
def unoptimized_mol():
    return a_molecule()
