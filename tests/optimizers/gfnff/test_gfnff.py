import numpy as np
import pytest
import stk

import stko
from tests.optimizers.utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)

pytest.importorskip("gfnff")


def _construction_bond_lengths(
    molecule: stk.ConstructedMolecule,
) -> list[float]:
    """Lengths of the bonds made during construction, in Å."""
    positions = molecule.get_position_matrix()
    return [
        float(
            np.linalg.norm(
                positions[info.get_bond().get_atom1().get_id()]
                - positions[info.get_bond().get_atom2().get_id()]
            )
        )
        for info in molecule.get_bond_infos()
        if info.get_building_block_id() is None
    ]


def test_gfnff_topo_pulls_in_construction_bonds(
    polymer: stk.ConstructedMolecule,
) -> None:
    # Built 3.15 Å long, too long for GFN-FF to find from the coordinates.
    assert min(_construction_bond_lengths(polymer)) > 3.0  # noqa: PLR2004
    optimizer = stko.OptimizerSequence(
        stko.GFNFFTopo(version="harmonic2020"),
        stko.GFNFFTopo(version="conformer2020"),
    )
    optimized = optimizer.optimize(polymer)

    is_equivalent_molecule(optimized, polymer)
    inequivalent_position_matrices(optimized, polymer)
    # C-C bonds between building blocks end near 1.53 Å.
    assert all(
        1.4 < length < 1.7  # noqa: PLR2004
        for length in _construction_bond_lengths(optimized)
    )


def test_gfnff_topo_without_bond_graph(
    polymer: stk.ConstructedMolecule,
) -> None:
    # GFN-FF's own bond perception doesn't see the long bonds.
    optimized = stko.GFNFFTopo(use_bond_graph=False).optimize(polymer)
    assert all(
        length > 3.0  # noqa: PLR2004
        for length in _construction_bond_lengths(optimized)
    )


def _trans_angles(molecule: stk.ConstructedMolecule) -> list[float]:
    """The two largest N-Pd-N angles at each Pd, in degrees."""
    positions = molecule.get_position_matrix()
    angles = []
    for atom in molecule.get_atoms():
        if atom.get_atomic_number() != 46:  # noqa: PLR2004
            continue
        centre = atom.get_id()
        ligands = [
            bond.get_atom2().get_id()
            if bond.get_atom1().get_id() == centre
            else bond.get_atom1().get_id()
            for bond in molecule.get_bonds()
            if centre in (bond.get_atom1().get_id(), bond.get_atom2().get_id())
        ]
        vectors = [positions[i] - positions[centre] for i in ligands]
        at_centre = sorted(
            np.degrees(
                np.arccos(
                    np.clip(
                        vectors[i]
                        @ vectors[j]
                        / np.linalg.norm(vectors[i])
                        / np.linalg.norm(vectors[j]),
                        -1,
                        1,
                    )
                )
            )
            for i in range(4)
            for j in range(i + 1, 4)
        )
        angles += at_centre[-2:]
    return angles


def _optimize(
    molecule: stk.ConstructedMolecule,
    *,
    square_planar: bool = False,
) -> stk.ConstructedMolecule:
    return stko.OptimizerSequence(
        stko.GFNFFTopo(
            version="harmonic2020", charge=4, square_planar=square_planar
        ),
        stko.GFNFFTopo(
            version="conformer2020", charge=4, square_planar=square_planar
        ),
    ).optimize(molecule)


def test_gfnff_topo_square_planar(lantern: stk.ConstructedMolecule) -> None:
    # Plain GFN-FF (the default) twists Pd(II) towards tetrahedral.
    free = _optimize(lantern)
    assert min(_trans_angles(free)) < 155  # noqa: PLR2004
    restrained = _optimize(lantern, square_planar=True)
    assert min(_trans_angles(restrained)) > 155  # noqa: PLR2004
