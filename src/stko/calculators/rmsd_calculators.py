import logging
import typing

import numpy as np
import stk
from scipy.spatial.distance import cdist

import stko
from stko.calculators.results.rmsd_results import RmsdResults
from stko.calculators.utilities import is_inequivalent_atom
from stko.utilities.exceptions import (
    DifferentAtomError,
    DifferentMoleculeError,
)

logger = logging.getLogger(__name__)


class RmsdCalculator:
    """
    Calculates the root mean square distance between molecules.

    This calculator will only work if the two molecules are the same
    and have the same atom ordering.

    No alignment of the two molecules occurs. However, both molecules
    are moved to a centroid position of (0, 0, 0).

    Examples:

        .. code-block:: python

            import stk
            import stko

            bb1 = stk.BuildingBlock('C1CCCCC1')
            calculator = stko.RmsdCalculator(bb1)
            results = calculator.get_results(stk.UFF().optimize(bb1))
            rmsd  = results.get_rmsd()

    """

    def __init__(
        self,
        initial_molecule: stk.Molecule,
        ignore_hydrogens: bool = False,
    ) -> None:
        """
        Initialize RMSD calculator with initial molecule.

        Parameters:

            initial_molecule:
                The :class:`.Molecule` to calculate RMSD from.

            ignore_hydrogens:
                ``True`` to ignore hydrogen atoms.

        """

        self._initial_molecule = initial_molecule
        self._ignore_hydrogens = ignore_hydrogens

    def _check_valid_comparison(self, mol: stk.Molecule) -> None:
        if mol.get_num_atoms() != (self._initial_molecule.get_num_atoms()):
            raise DifferentMoleculeError(
                f"{self._initial_molecule} and {mol} are not "
                "equivalent with different numbers of atoms."
            )

        smiles1 = stk.Smiles().get_key(self._initial_molecule)
        smiles2 = stk.Smiles().get_key(mol)
        if smiles1 != smiles2:
            raise DifferentMoleculeError(
                f"{self._initial_molecule} and {mol} are not "
                "equivalent with different smiles strings."
            )

        atoms1 = self._initial_molecule.get_atoms()
        atoms2 = mol.get_atoms()
        for atom1, atom2 in zip(atoms1, atoms2):
            if is_inequivalent_atom(atom1, atom2):
                raise DifferentAtomError(
                    f"{atom1} and {atom2} are not equivalent."
                )

    def _calculate_rmsd(self, mol: stk.Molecule) -> float:
        if self._ignore_hydrogens:
            initial_atom_ids = (
                i.get_id()
                for i in self._initial_molecule.get_atoms()
                if i.get_atomic_number() != 1
            )
            mol_atom_ids = (
                i.get_id()
                for i in mol.get_atoms()
                if i.get_atomic_number() != 1
            )
        else:
            initial_atom_ids = None
            mol_atom_ids = None

        initial_atom_positions = self._initial_molecule.get_atomic_positions(
            atom_ids=initial_atom_ids,
        )
        mol_atom_positions = mol.get_atomic_positions(
            atom_ids=mol_atom_ids,
        )

        pos_mat1 = np.array(list(initial_atom_positions))
        pos_mat2 = np.array(list(mol_atom_positions))

        deviations = pos_mat1 - pos_mat2
        N = len(pos_mat1)
        return np.sqrt(np.sum(deviations * deviations) / N)

    def calculate(self, mol: stk.Molecule) -> typing.Iterable[float]:
        self._check_valid_comparison(mol)
        self._initial_molecule = self._initial_molecule.with_centroid(
            position=np.array((0, 0, 0)),
        )
        mol = mol.with_centroid(np.array((0, 0, 0)))
        yield self._calculate_rmsd(mol)

    def get_results(self, mol: stk.Molecule) -> stko.RmsdResults:
        """
        Calculate the RMSD between `mol` and the initial molecule.

        Parameters:

            mol:
                The :class:`.Molecule` to calculate RMSD to.

        Returns:

            The RMSD between the molecules.

        """

        return RmsdResults(self.calculate(mol))


class RmsdMappedCalculator(RmsdCalculator):
    """
    Calculates the root mean square distance between molecules.

    This calculator allows for different molecules but they should be
    aligned, see the example below. It will calculate the RMSD based on
    the nearest atom of the same class (element). Both molecules are
    moved to a centroid position of (0, 0, 0). The number of atoms is
    based on the `mol` input into `calculate`.

    Warning: the RMSD depends on the order, i.e. it is not guaranteed
    to be the same when you switch the initial and test molecule.

    Examples:

        .. code-block:: python

            import stk
            import stko
            import numpy as np

            bb1 = stk.BuildingBlock('C1CCCCC1')
            # Fake rotation of new molecule.
            bb2 = stk.BuildingBlock('C1CCCCC1').with_rotation_about_axis(
                1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
            )

            # Get RMSD without alignment.
            calculator = stko.RmsdMappedCalculator(bb1)
            results = calculator.get_results(bb2)
            rmsd  = results.get_rmsd()

            # Align the molecules.
            optimizer = stko.Aligner(bb1, (('C', 'C'), ))
            aligned_bb2 = optimizer.optimize(bb2)

            calculator = stko.RmsdMappedCalculator(bb1)
            results = calculator.get_results(aligned_bb2)
            rmsd  = results.get_rmsd()

    """

    def _calculate_rmsd(self, mol: stk.Molecule) -> float:
        if self._ignore_hydrogens:
            initial_atom_ids = (
                i.get_id()
                for i in self._initial_molecule.get_atoms()
                if i.get_atomic_number() != 1
            )
            mol_atom_ids = (
                i.get_id()
                for i in mol.get_atoms()
                if i.get_atomic_number() != 1
            )
        else:
            initial_atom_ids = None
            mol_atom_ids = None

        initial_atom_positions = self._initial_molecule.get_atomic_positions(
            atom_ids=initial_atom_ids,
        )
        initial_atoms = tuple(
            self._initial_molecule.get_atoms(atom_ids=initial_atom_ids)
        )
        mol_atoms = tuple(mol.get_atoms(atom_ids=mol_atom_ids))
        atom_matrix = np.zeros((len(initial_atoms), len(mol_atoms)))
        for i in range(len(mol_atoms)):
            a1_num = mol_atoms[i].get_atomic_number()
            for j in range(len(initial_atoms)):
                a2_num = initial_atoms[j].get_atomic_number()
                atom_matrix[j][i] = a1_num == a2_num

        mol_atom_positions = mol.get_atomic_positions(
            atom_ids=mol_atom_ids,
        )
        pos_mat1 = np.array(list(initial_atom_positions))
        pos_mat2 = np.array(list(mol_atom_positions))

        N = len(pos_mat2)
        distances = cdist(pos_mat1, pos_mat2)
        new_array = np.where(atom_matrix, distances, 1e24)
        # Handle situation where one molecule does not have an atom
        # type that the other does.
        deviations = np.array(
            [i for i in np.amin(new_array, axis=1) if i != 1e24]
        )
        return np.sqrt(np.sum(deviations * deviations) / N)

    def calculate(self, mol: stk.Molecule) -> typing.Iterable[float]:
        self._initial_molecule = self._initial_molecule.with_centroid(
            position=np.array((0, 0, 0)),
        )
        mol = mol.with_centroid(np.array((0, 0, 0)))
        yield self._calculate_rmsd(mol)
