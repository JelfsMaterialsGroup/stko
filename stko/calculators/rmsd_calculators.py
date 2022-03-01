"""
RMSD Calculators
=================

#. :class:`.RmsdCalculator`

Calculator of Root Mean Square Distance between two molecules.

"""

import logging
import stk
import numpy as np

from .calculators import Calculator
from .results import RmsdResults
from ..utilities import is_inequivalent_atom


logger = logging.getLogger(__name__)


class RmsdCalculatorError(Exception):
    ...


class DifferentAtomException(RmsdCalculatorError):
    ...


class DifferentMoleculeException(RmsdCalculatorError):
    ...


class RmsdCalculator(Calculator):
    """
    Calculates the root mean square distance between molecules.

    This calculator will only work if the two molecules are the same
    and have the same atom ordering.

    No alignment of the two molecules occurs. However, both molecules
    are moved to a centroid position of (0, 0, 0).

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        bb1 = stk.BuildingBlock('C1CCCCC1')
        calculator = stko.RmsdCalculator(bb1)
        results = calculator.get_results(stk.UFF().optimize(bb1))
        rmsd  = results.get_rmsd()

    """

    def __init__(self, initial_molecule, ignore_hydrogens=False):
        """
        Initialize RMSD calculator with initial molecule.

        Parameters
        ----------
        initial_molecule : :class:`.Molecule`
            The :class:`.Molecule` to calculate RMSD from.

        ignore_hydrogens : :class:`bool`, optional
            ``True`` to ignore hydrogen atoms.

        """

        self._initial_molecule = initial_molecule
        self._ignore_hydrogens = ignore_hydrogens

    def _check_valid_comparison(self, mol):

        if mol.get_num_atoms() != (
            self._initial_molecule.get_num_atoms()
        ):
            raise DifferentMoleculeException(
                f'{self._initial_molecule} and {mol} are not '
                'equivalent with different numbers of atoms.'
            )

        smiles1 = stk.Smiles().get_key(self._initial_molecule)
        smiles2 = stk.Smiles().get_key(mol)
        if smiles1 != smiles2:
            raise DifferentMoleculeException(
                f'{self._initial_molecule} and {mol} are not '
                'equivalent with different smiles strings.'
            )

        atoms1 = self._initial_molecule.get_atoms()
        atoms2 = mol.get_atoms()
        for atom1, atom2 in zip(atoms1, atoms2):
            if is_inequivalent_atom(atom1, atom2):
                raise DifferentAtomException(
                    f'{atom1} and {atom2} are not equivalent.'
                )

    def _calculate_rmsd(self, mol):
        if self._ignore_hydrogens:
            initial_atom_ids = (
                i.get_id() for i in self._initial_molecule.get_atoms()
                if i.get_atomic_number() != 1
            )
            mol_atom_ids = (
                i.get_id() for i in mol.get_atoms()
                if i.get_atomic_number() != 1
            )
        else:
            initial_atom_ids = None
            mol_atom_ids = None

        initial_atom_positions = (
            self._initial_molecule.get_atomic_positions(
                atom_ids=initial_atom_ids,
            )
        )
        mol_atom_positions = (
            mol.get_atomic_positions(
                atom_ids=mol_atom_ids,
            )
        )

        pos_mat1 = np.array(list(initial_atom_positions))
        pos_mat2 = np.array(list(mol_atom_positions))

        deviations = pos_mat1 - pos_mat2
        N = len(pos_mat1)
        return np.sqrt(np.sum(deviations * deviations) / N)

    def calculate(self, mol):
        self._check_valid_comparison(mol)
        self._initial_molecule = self._initial_molecule.with_centroid(
            position=(0, 0, 0)
        )
        mol = mol.with_centroid((0, 0, 0))
        yield self._calculate_rmsd(mol)

    def get_results(self, mol):
        """
        Calculate the RMSD between `mol` and the initial molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` to calculate RMSD to.

        Returns
        -------
        :class:`.RmsdResults`
            The RMSD between the molecules.

        """

        return RmsdResults(self.calculate(mol))
