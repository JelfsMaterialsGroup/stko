import pytest
import sys

from stko import GulpUFFOptimizer, GulpUFFMDOptimizer


gulp = pytest.mark.skipif(
    all('gulp_path' not in x for x in sys.argv),
    reason="Only run when explicitly asked."
)


@gulp
def test_optimizer1(gulp_path):
    gulpuffoptimizer = GulpUFFOptimizer(
        gulp_path=gulp_path,
        maxcyc=1000,
        metal_FF=None,
        metal_ligand_bond_order=None,
        conjugate_gradient=False,
        periodic=False,
        output_dir=None,
    )


    assert True


@gulp
def test_optimizer2(gulp_path):
    gulpuffmdoptimizer = GulpUFFMDOptimizer(
        gulp_path=gulp_path,
        metal_FF=None,
        metal_ligand_bond_order=None,
        output_dir=None,
        integrator='stochastic',
        ensemble='nvt',
        temperature=300,
        equilbration=1.0,
        production=10.0,
        timestep=1.0,
        N_conformers=10,
        opt_conformers=True,
        save_conformers=False,
    )

    assert True
