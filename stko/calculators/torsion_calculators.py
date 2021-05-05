"""
Torsion Calculators
===================

#. :class:`.TorsionCalculator`

"""

import logging

from .calculators import Calculator
from .results import TorsionResults, ConstructedMoleculeTorsionResults
from rdkit.Chem import TorsionFingerprints
from ..molecular.torsion import Torsion

logger = logging.getLogger(__name__)


class TorsionCalculator(Calculator):
    """
    Uses rdkit to extract all torsions in a molecule.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        # Create a molecule whose energy we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the calculator.
        tc = stko.TorsionCalculator()

        # Extract the torsions.
        ........


    """

    def calculate(self, mol):
        yield tuple(
            Torsion(*mol.get_atoms(atoms[0]))
            for atoms, _ in (
                TorsionFingerprints.CalculateTorsionLists(
                    mol.to_rdkit_mol()
                )[0]
            )
        )

    def get_results(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`.TorsionResults`
            The torsions of the molecule.

        """

        return TorsionResults(self.calculate(mol), mol)


class ConstructedMoleculeTorsionCalculator(TorsionCalculator):
    """
    Uses rdkit to extract all torsions in a molecule.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        # Create a molecule whose energy we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the calculator.
        mmff = stko.TorsionCalculator()

        # Extract the torsions.
        ........


    """

    def get_results(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`.TorsionResults`
            The torsions of the molecule.

        """

        return ConstructedMoleculeTorsionResults(
            generator=self.calculate(mol),
            mol=mol,
        )
