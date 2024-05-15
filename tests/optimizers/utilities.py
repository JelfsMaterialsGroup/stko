import itertools as it

import numpy as np
import stk


def is_equivalent_atom(atom1: stk.Atom, atom2: stk.Atom) -> None:
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
    assert atom1.get_atomic_number() == atom2.get_atomic_number()
    assert atom1.__class__ is atom2.__class__


def is_equivalent_bond(bond1: stk.Bond, bond2: stk.Bond) -> None:
    assert bond1.__class__ is bond2.__class__
    assert bond1.get_order() == bond2.get_order()
    assert bond1.get_periodicity() == bond2.get_periodicity()
    is_equivalent_atom(bond1.get_atom1(), bond2.get_atom1())
    is_equivalent_atom(bond1.get_atom2(), bond2.get_atom2())


def is_equivalent_molecule(
    molecule1: stk.Molecule, molecule2: stk.Molecule
) -> None:
    atoms = it.zip_longest(
        molecule1.get_atoms(),
        molecule2.get_atoms(),
    )
    for atom1, atom2 in atoms:
        is_equivalent_atom(atom1, atom2)

    bonds = it.zip_longest(
        molecule1.get_bonds(),
        molecule2.get_bonds(),
    )
    for bond1, bond2 in bonds:
        is_equivalent_bond(bond1, bond2)


def inequivalent_position_matrices(
    molecule1: stk.Molecule, molecule2: stk.Molecule
) -> None:
    pos1 = molecule1.get_position_matrix()
    pos2 = molecule2.get_position_matrix()

    assert not np.allclose(pos1, pos2)
