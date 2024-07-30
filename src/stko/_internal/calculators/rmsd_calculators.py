import logging

import numpy as np
import stk
from rmsd import kabsch_rmsd
from scipy.spatial.distance import cdist

from stko._internal.calculators.results.rmsd_results import RmsdResults
from stko._internal.calculators.utilities import is_inequivalent_atom
from stko._internal.utilities.exceptions import (
    DifferentAtomError,
    DifferentMoleculeError,
)

logger = logging.getLogger(__name__)


class RmsdCalculator:
    """Calculates the root mean square distance between molecules.

    This calculator will only work if the two molecules are the same
    and have the same atom ordering.

    No alignment of the two molecules occurs. However, both molecules
    are moved to a centroid position of (0, 0, 0).

    Parameters:
        initial_molecule:
            The :class:`stk.Molecule` to calculate RMSD from.

        ignore_hydrogens:
            ``True`` to ignore hydrogen atoms.


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
        self._initial_molecule = initial_molecule
        self._ignore_hydrogens = ignore_hydrogens

    def _check_valid_comparison(self, mol: stk.Molecule) -> None:
        if mol.get_num_atoms() != (self._initial_molecule.get_num_atoms()):
            msg = (
                f"{self._initial_molecule} and {mol} are not "
                "equivalent with different numbers of atoms."
            )
            raise DifferentMoleculeError(msg)

        smiles1 = stk.Smiles().get_key(self._initial_molecule)
        smiles2 = stk.Smiles().get_key(mol)
        if smiles1 != smiles2:
            msg = (
                f"{self._initial_molecule} and {mol} are not "
                "equivalent with different smiles strings."
            )
            raise DifferentMoleculeError(msg)

        atoms1 = self._initial_molecule.get_atoms()
        atoms2 = mol.get_atoms()
        for atom1, atom2 in zip(atoms1, atoms2, strict=False):
            if is_inequivalent_atom(atom1, atom2):
                msg = f"{atom1} and {atom2} are not equivalent."
                raise DifferentAtomError(msg)

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
        n = len(pos_mat1)
        return np.sqrt(np.sum(deviations * deviations) / n)

    def calculate(self, mol: stk.Molecule) -> float:
        self._check_valid_comparison(mol)
        self._initial_molecule = self._initial_molecule.with_centroid(
            position=np.array((0, 0, 0)),
        )
        mol = mol.with_centroid(np.array((0, 0, 0)))
        return self._calculate_rmsd(mol)

    def get_results(self, mol: stk.Molecule) -> RmsdResults:
        """Calculate the RMSD between `mol` and the initial molecule.

        Parameters:
            mol:
                The :class:`stk.Molecule` to calculate RMSD to.

        Returns:
            The RMSD between the molecules.

        """
        return RmsdResults(self.calculate(mol))


class RmsdMappedCalculator(RmsdCalculator):
    """Calculates the root mean square distance between molecules.

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

        n = len(pos_mat2)
        distances = cdist(pos_mat1, pos_mat2)
        initial_value = 1e24
        new_array = np.where(atom_matrix, distances, initial_value)
        # Handle situation where one molecule does not have an atom
        # type that the other does.
        deviations = np.array(
            [i for i in np.amin(new_array, axis=1) if i != initial_value]
        )
        return np.sqrt(np.sum(deviations * deviations) / n)

    def calculate(self, mol: stk.Molecule) -> float:
        self._initial_molecule = self._initial_molecule.with_centroid(
            position=np.array((0, 0, 0)),
        )
        mol = mol.with_centroid(np.array((0, 0, 0)))
        return self._calculate_rmsd(mol)


class KabschRmsdCalculator:
    """Calculates the root mean square distance between molecules.

    This calculator uses the rmsd package with the default settings and no
    reordering.

    See Also:
        * rmsd https://github.com/charnley/rmsd

    This calculator will only work if the two molecules are the same
    and have the same atom ordering.

    Parameters:
        initial_molecule:
            The :class:`stk.Molecule` to calculate RMSD from.

    Examples:
        .. code-block:: python

            import stk
            import stko

            bb1 = stk.BuildingBlock('C1CCCCC1')
            calculator = stko.KabschRmsdCalculator(bb1)
            results = calculator.get_results(stk.UFF().optimize(bb1))
            rmsd  = results.get_rmsd()

    """

    def __init__(self, initial_molecule: stk.Molecule) -> None:
        self._initial_molecule = initial_molecule

    def _calculate_rmsd(self, mol: stk.Molecule) -> float:
        p_coord = self._initial_molecule.get_position_matrix()
        q_coord = mol.get_position_matrix()
        return kabsch_rmsd(p_coord, q_coord)

    def calculate(self, mol: stk.Molecule) -> float:
        self._initial_molecule = self._initial_molecule.with_centroid(
            position=np.array((0, 0, 0)),
        )
        mol = mol.with_centroid(np.array((0, 0, 0)))
        return self._calculate_rmsd(mol)

    def get_results(self, mol: stk.Molecule) -> RmsdResults:
        """Calculate the RMSD between `mol` and the initial molecule.

        Parameters:
            mol:
                The :class:`stk.Molecule` to calculate RMSD to.

        Returns:
            The RMSD between the molecules.

        """
        return RmsdResults(self.calculate(mol))
