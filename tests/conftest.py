

def pytest_addoption(parser):
    parser.addoption('--macromodel_path', default='')


def pytest_generate_tests(metafunc):
    if 'macromodel_path' in metafunc.fixturenames:
        macromodel_path = metafunc.config.getoption('macromodel_path')
        metafunc.parametrize('macromodel_path', [macromodel_path])
