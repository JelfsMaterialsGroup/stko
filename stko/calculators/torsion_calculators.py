"""
Torsion Calculators
===================

#. :class:`.TorsionCalculator`
#. :class:`.ConstructedMoleculeTorsionCalculator`

Methods to extract torsions from a molecule or constructed molecule.

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

    Note that the rdkit [1]_ function we use only outputs
    one torsion for each rotatable bond. We use the
    `TorsionFingerprints.CalculateTorsionLists` method.

    Examples
    --------

    .. code-block:: python

        import stk
        import stko

        # Create a molecule whose torsions we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the calculator.
        tc = stko.TorsionCalculator()

        # Extract the torsions.
        tc_results = tc.get_results(mol1)
        for t, ang in tc_results.get_torsion_angles():
            print(t, ang, t.get_atom_ids())

    References
    ----------
    .. [1] http://rdkit.org/docs/source/
    rdkit.Chem.TorsionFingerprints.html

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

    Note that the rdkit [1]_ function we use only outputs
    one torsion for each rotatable bond. We use the
    `TorsionFingerprints.CalculateTorsionLists` method.

    Examples
    --------

    .. code-block:: python

        import stk
        import stko

        # Create a molecule whose energy we want to know.
        bb1 = stk.BuildingBlock('NCCNCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
        polymer = stk.ConstructedMolecule(
            stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit="AB",
                orientations=[0, 0],
                num_repeating_units=1
            )
        )

        # Create the calculator.
        tc = stko.ConstructedMoleculeTorsionCalculator()

        # Extract the torsions.
        tc_results = tc.get_results(polymer)

        # Get information about torsions in building blocks and in the
        # ConstructedMolecule.
        for t in tc_results.get_torsion_infos():
            print(
                'c', t.get_torsion(),
                t.get_building_block(),
                t.get_building_block_id(),
                t.get_building_block_torsion(),
            )

    References
    ----------
    .. [1] http://rdkit.org/docs/source/
    rdkit.Chem.TorsionFingerprints.html

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
