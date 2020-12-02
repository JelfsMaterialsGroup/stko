import stk
import stko
import pytest
import os
import stk

@pytest.fixture
def benzene_build():
    """
    Benzene fixture with distorted geometry
    """
    path_to_current_file = os.path.realpath(__file__)
    current_directory = os.path.split(path_to_current_file)[0]
    benzene = os.path.join(current_directory, "benzene.mol")
    mol = stk.BuildingBlock.init_from_file(benzene)
    return mol