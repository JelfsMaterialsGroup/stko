import pytest
import os
import stk


@pytest.fixture
def benzene_build():
    """
    Benzene fixture with distorted geometry.

    """

    path_to_current_file = os.path.realpath(__file__)
    current_directory = os.path.split(path_to_current_file)[0]
    benzene = os.path.join(current_directory, "benzene.mol")
    mol = stk.BuildingBlock.init_from_file(benzene)

    return mol


def pytest_addoption(parser):
    parser.addoption('--gulp_path', default='')
    parser.addoption('--macromodel_path', default='')


def pytest_generate_tests(metafunc):
    if 'gulp_path' in metafunc.fixturenames:
        gulp_path = metafunc.config.getoption('gulp_path')
        metafunc.parametrize('gulp_path', [gulp_path])
    elif 'macromodel_path' in metafunc.fixturenames:
        macromodel_path = metafunc.config.getoption('macromodel_path')
        metafunc.parametrize('macromodel_path', [macromodel_path])
