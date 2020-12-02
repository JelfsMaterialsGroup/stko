def pytest_addoption(parser):
    parser.addoption('--gulp_path', default='')


def pytest_generate_tests(metafunc):
    if 'gulp_path' in metafunc.fixturenames:
        gulp_path = metafunc.config.getoption('gulp_path')
        metafunc.parametrize('gulp_path', [gulp_path])
