import pytest
import sys
import os
from os.path import join

from stko import GulpUFFOptimizer, GulpUFFMDOptimizer
from.utilities import compare_benzenes


odir = join(
    os.path.dirname(os.path.abspath(__file__)),
    'gulp_tests_output',
)
if not os.path.exists(odir):
    os.mkdir(odir)

gulp = pytest.mark.skipif(
    all('gulp_path' not in x for x in sys.argv),
    reason="Only run when explicitly asked."
)


@gulp
def test_optimizer1(gulp_path, benzene_build):
    gulpuffoptimizer = GulpUFFOptimizer(
        gulp_path=gulp_path,
        maxcyc=1000,
        metal_FF=None,
        metal_ligand_bond_order=None,
        conjugate_gradient=False,
        periodic=False,
        output_dir=join(odir, 'test_optimizer1'),
    )
    gulpuffoptimizer.assign_FF(benzene_build)
    opt_benzene = gulpuffoptimizer.optimize(benzene_build)
    compare_benzenes(
        initial_molecule=benzene_build,
        optimized_molecule=opt_benzene,
    )


@gulp
def test_optimizer2(gulp_path, benzene_build):
    gulpuffmdoptimizer = GulpUFFMDOptimizer(
        gulp_path=gulp_path,
        metal_FF=None,
        metal_ligand_bond_order=None,
        output_dir=join(odir, 'test_optimizer2'),
        integrator='stochastic',
        ensemble='nvt',
        temperature=300,
        equilbration=1.0,
        production=5.0,
        timestep=0.5,
        N_conformers=4,
        opt_conformers=True,
        save_conformers=False,
    )
    gulpuffmdoptimizer.assign_FF(benzene_build)
    opt_benzene = gulpuffmdoptimizer.optimize(benzene_build)
    compare_benzenes(
        initial_molecule=benzene_build,
        optimized_molecule=opt_benzene,
    )


@gulp
def test_optimizer3(gulp_path, benzene_build):
    gulpuffmdoptimizer = GulpUFFMDOptimizer(
        gulp_path=gulp_path,
        metal_FF=None,
        metal_ligand_bond_order=None,
        output_dir=join(odir, 'test_optimizer3'),
        integrator='stochastic',
        ensemble='nvt',
        temperature=150.0,
        equilbration=1.0,
        production=10.0,
        timestep=0.5,
        N_conformers=4,
        opt_conformers=False,
        save_conformers=False,
    )
    gulpuffmdoptimizer.assign_FF(benzene_build)
    opt_benzene = gulpuffmdoptimizer.optimize(benzene_build)
    compare_benzenes(
        initial_molecule=benzene_build,
        optimized_molecule=opt_benzene,
    )
