import pytest
import numpy as np
import stko
import stk


def a_molecule():
    return stk.BuildingBlock(smiles='CCCCCC')


@pytest.fixture
def unoptimized_mol():
    return a_molecule()
