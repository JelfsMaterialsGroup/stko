import pytest
import numpy as np
from stko import GulpUFFOptimizer, GulpUFFMDOptimizer, get_metal_atoms
from .conftest import a_molecule


class FakeGulpUFFOptimizer(GulpUFFOptimizer):

    def optimize(self, mol):
        return a_molecule().with_centroid(np.array(([1, 3, 3])))


class FakeGulpUFFMDOptimizer(GulpUFFMDOptimizer):

    def optimize(self, mol):
        return a_molecule().with_centroid(np.array(([1, 3, 3])))


@pytest.fixture
def position_section():
    return (
        '\ncartesian\n'
        'C1 core -0.74031 0.03171 0.10194\n'
        'C1 core 0.75974 -0.01804 -0.03434\n'
        'H1 core -1.14779 -0.63142 -0.70939\n'
        'H1 core -1.11274 -0.39237 1.04655\n'
        'H1 core -1.12625 1.0535 -0.10345\n'
        'H1 core 0.97939 0.03507 -1.13219\n'
        'H1 core 1.25573 0.87341 0.40438\n'
        'H1 core 1.13222 -0.95187 0.4265\n'
    )


@pytest.fixture
def bond_section():
    return (
        '\nconnect 1 2 \n'
        'connect 1 3 \n'
        'connect 1 4 \n'
        'connect 1 5 \n'
        'connect 2 6 \n'
        'connect 2 7 \n'
        'connect 2 8 \n'
    )


@pytest.fixture
def species_section():
    return '\nspecies\nC1 C_3\nH1 H_\n'


def test_gulp_position_section(unoptimized_mol, position_section):
    opt = FakeGulpUFFOptimizer(
        gulp_path='',
        maxcyc=1000,
        metal_FF=None,
        metal_ligand_bond_order=None,
        conjugate_gradient=False,
        output_dir='',
    )
    opt.assign_FF(unoptimized_mol)
    type_translator = opt._type_translator()
    test = opt._position_section(unoptimized_mol, type_translator)
    assert position_section == test


def test_gulp_bond_section(unoptimized_mol, bond_section):
    opt = FakeGulpUFFOptimizer(
        gulp_path='',
        maxcyc=1000,
        metal_FF=None,
        metal_ligand_bond_order=None,
        conjugate_gradient=False,
        output_dir='',
    )
    opt.assign_FF(unoptimized_mol)
    metal_atoms = get_metal_atoms(unoptimized_mol)
    test = opt._bond_section(unoptimized_mol, metal_atoms)
    assert bond_section == test


def test_gulp_species_section(unoptimized_mol, species_section):
    opt = FakeGulpUFFOptimizer(
        gulp_path='',
        maxcyc=1000,
        metal_FF=None,
        metal_ligand_bond_order=None,
        conjugate_gradient=False,
        output_dir='',
    )
    opt.assign_FF(unoptimized_mol)
    type_translator = opt._type_translator()
    test = opt._species_section(type_translator)
    assert species_section == test
