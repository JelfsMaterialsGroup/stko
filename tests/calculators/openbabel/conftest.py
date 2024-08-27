import pytest
import stk


def a_molecule() -> stk.BuildingBlock:
    return stk.BuildingBlock(smiles="CCCCCC")


@pytest.fixture
def unoptimized_mol() -> stk.BuildingBlock:
    return a_molecule()
