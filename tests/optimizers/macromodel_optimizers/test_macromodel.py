import pytest
import sys
from stko import MacroModelMD, MacroModelForceField


macromodel_check = pytest.mark.skipif(
    all('macromodel_path' not in x for x in sys.argv),
    reason="Only run when MacroModel path is explicitly given."
)


@pytest.fixture(scope="function")
def macromodel_md(macromodel_path):
    return MacroModelMD(macromodel_path)


@pytest.fixture(scope="function")
def macromodel_forcefield(tmpdir, macromodel_path):
    return MacroModelForceField(
        force_field=16,
        output_dir=tmpdir,
        macromodel_path=macromodel_path
    )


@macromodel_check
class TestMacroModelFF:
    def test_generate_com(
        self,
        macromodel_forcefield,
        tmp_path,
        benzene_build
    ):
        # Write the ``.com`` file.
        run_name = 'test_com'
        macromodel_forcefield._generate_com(
            benzene_build, str(tmp_path/run_name)
        )
        with open(str(tmp_path / (run_name+'.com')), 'r') as f:
            com = f.readlines()
        com_path = (
            '/home/sbennett/PhD/stko/tests/optimizers/'
            'macromodel_optimizers/benzene_test.com'
        )
        with open(com_path, 'r') as f:
            test_com = f.readlines()[2:]

        # Do not compare first two lines, as these correspond
        assert test_com == com[2:]
