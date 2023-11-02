"""
Aligner
=======

Optimizer that aligns two molecules using `Spindry`.

"""

import logging
from itertools import product

import stk
import numpy as np
import spindry as spd
from scipy.spatial.distance import cdist

from stko.calculators.rmsd_calculators import RmsdMappedCalculator
from stko.optimizers.optimizers import Optimizer

logger = logging.getLogger(__name__)


class AlignmentPotential(spd.Potential):
    def __init__(self, matching_pairs: tuple[tuple[str]], width: float):
        self._matching_pairs = matching_pairs
        self._width = width

    def _potential(self, distance: np.ndarray) -> np.ndarray:
        return self._width * (distance**2)

    def _combine_atoms(
        self,
        atoms1: tuple[spd.Atom, ...],
        atoms2: tuple[spd.Atom, ...],
    ) -> np.ndarray:
        len1 = len(atoms1)
        len2 = len(atoms2)

        mixed = np.zeros((len1, len2))
        for i in range(len1):
            for j in range(len2):
                a1e = atoms1[i].get_element_string()
                a2e = atoms2[j].get_element_string()
                if tuple(sorted((a1e, a2e))) in self._matching_pairs:
                    mixed[i, j] = True
                else:
                    mixed[i, j] = False

        return mixed

    def compute_potential(self, supramolecule: spd.SupraMolecule) -> float:
        component_position_matrices = list(
            i.get_position_matrix() for i in supramolecule.get_components()
        )
        component_atoms = list(
            tuple(j for j in i.get_atoms())
            for i in supramolecule.get_components()
        )
        pair_dists = cdist(
            component_position_matrices[0],
            component_position_matrices[1],
        )
        sigmas = self._combine_atoms(component_atoms[0], component_atoms[1])

        new_array = np.where(sigmas, pair_dists, 1e24)
        distances = np.array(
            [i for i in np.amin(new_array, axis=1) if i != 1e24]
        )
        return np.sum(self._potential(distance=distances))


class Aligner(Optimizer):
    """
    Use SpinDry to align two molecules.[1]_

    Examples:

        .. code-block:: python

            import stk
            import stko
            import numpy as np


            # For this example, we have produced rotated molecules.
            mol = stk.BuildingBlock('NCCNCCN')
            mol2 = mol.with_rotation_about_axis(
                1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
            )
            aligner = stko.Aligner(mol2, (('N', 'N'), ))
            mol = aligner.optimize(mol)

    References:

        .. [1] https://github.com/andrewtarzia/SpinDry

    """

    def __init__(
        self,
        initial_molecule: stk.Molecule,
        matching_pairs: tuple[tuple[str]],
    ) -> None:
        """
        Initialize aligner optimizer.

        Parameters:

            initial_molecule:
                Molecule to align to.

            matching_pairs:
                Pairs of atom types to use in alignment.

        """

        self._initial_molecule = initial_molecule.with_centroid(
            np.array((0, 0, 0)),
        )
        self._matching_pairs = matching_pairs

    def _get_supramolecule(self, mol: stk.Molecule) -> spd.SupraMolecule:
        host_molecule = spd.Molecule(
            atoms=(
                spd.Atom(
                    id=atom.get_id(),
                    element_string=atom.__class__.__name__,
                )
                for atom in self._initial_molecule.get_atoms()
            ),
            bonds=(
                spd.Bond(
                    id=i,
                    atom_ids=(
                        bond.get_atom1().get_id(),
                        bond.get_atom2().get_id(),
                    ),
                )
                for i, bond in enumerate(self._initial_molecule.get_bonds())
            ),
            position_matrix=(self._initial_molecule.get_position_matrix()),
        )
        guest_molecule = spd.Molecule(
            atoms=(
                spd.Atom(
                    id=atom.get_id(),
                    element_string=atom.__class__.__name__,
                )
                for atom in mol.get_atoms()
            ),
            bonds=(
                spd.Bond(
                    id=i,
                    atom_ids=(
                        bond.get_atom1().get_id(),
                        bond.get_atom2().get_id(),
                    ),
                )
                for i, bond in enumerate(mol.get_bonds())
            ),
            position_matrix=mol.get_position_matrix(),
        )

        return spd.SupraMolecule.init_from_components(
            components=(host_molecule, guest_molecule),
        )

    def _align_molecules(self, mol: stk.Molecule) -> stk.Molecule:
        supramolecule = self._get_supramolecule(mol)
        for comp in supramolecule.get_components():
            pass
        mol = mol.with_position_matrix(comp.get_position_matrix())
        cg = spd.Spinner(
            step_size=0.2,
            rotation_step_size=0.2,
            num_conformers=50,
            max_attempts=250,
            potential_function=AlignmentPotential(
                matching_pairs=self._matching_pairs,
                width=2,
            ),
        )
        conformer = cg.get_final_conformer(
            supramolecule=supramolecule,
            movable_components=(1,),
        )

        for comp in conformer.get_components():
            pass
        mol = mol.with_position_matrix(comp.get_position_matrix())

        return mol

    def optimize(self, mol: stk.Molecule) -> stk.Molecule:
        rotation_axes = (
            None,
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            (1, 0, 1),
            (1, 1, 0),
            (0, 1, 1),
        )
        angles = (
            np.radians(60),
            np.radians(90),
            np.radians(120),
            np.radians(180),
            np.radians(240),
            np.radians(270),
        )

        min_rmsd = 1e10
        for r, rot in enumerate(product(rotation_axes, angles)):
            aligned_mol = mol.with_centroid(np.array((0, 0, 0)))
            if rot[0] is not None:
                aligned_mol = aligned_mol.with_rotation_about_axis(
                    angle=rot[1],
                    axis=np.array(rot[0]),
                    origin=np.array((0, 0, 0)),
                )
            elif r != 0:
                continue

            aligned_mol = aligned_mol.with_centroid(np.array((0, 0, 0)))
            aligned_mol = self._align_molecules(aligned_mol)
            rmsd_calculator = RmsdMappedCalculator(self._initial_molecule)
            rmsd = rmsd_calculator.get_results(aligned_mol).get_rmsd()
            if rmsd < min_rmsd:
                min_rmsd = rmsd
                final_mol = aligned_mol.with_centroid(np.array((0, 0, 0)))

        return final_mol
