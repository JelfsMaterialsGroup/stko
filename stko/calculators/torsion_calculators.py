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
from ..utilities.utilities import get_atom_maps

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
        Calculate the torsions of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose torsions are to be calculated.

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
        Calculate the torsions of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose torsions are to be calculated.

        Returns
        -------
        :class:`.TorsionResults`
            The torsions of the molecule.

        """

        return ConstructedMoleculeTorsionResults(
            generator=self.calculate(mol),
            mol=mol,
        )


class MatchedTorsionCalculator(ConstructedMoleculeTorsionCalculator):
    """
    Matches rdkit generated torsions with building block torsions.

    """

    def calculate(self, mol):
        """
        Extract torsions with rdkit, then match to building blocks.

        This method loops through each rdkit generated torsion. For
        each torsion, it checks if the two interior atoms of the
        torsion come from the two interior (central) atoms of some
        torsion in a building block. If so, it replaces the end atoms
        of the torsion with the atoms corresponding to the end atoms
        of the underlying building block torsion.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` whose torsions are to be
            calculated.

        """

        torsions = list(next(super().calculate(mol)))
        atom_maps = get_atom_maps(mol)

        # Loop over torsions, updating each to match a building block
        # torsion if possible.
        for i, torsion in enumerate(torsions):
            # Atom ids, atom infos, and building block ids of current
            # torsion in constructed molecule.
            atom_ids = list(torsion.get_atom_ids())
            atom_infos = list(mol.get_atom_infos(atom_ids))
            build_block_ids = [
                atom_info.get_building_block_id()
                for atom_info in atom_infos
            ]

            if build_block_ids[1] is None:
                continue
            # Check if two central atoms of torsion are from the same
            # building block. If not, leave this torsion alone.
            if (build_block_ids[1] != build_block_ids[2]):
                continue

            build_block_torsions = TorsionCalculator().get_results(
                atom_infos[1].get_building_block()
            ).get_torsions()
            atom_map = atom_maps[build_block_ids[1]]

            # Look for a torsion in the building block that has the
            # same central atoms.
            for bb_torsion in build_block_torsions:
                try:
                    matched_atoms = [
                        atom_map[atom_id]
                        for atom_id in bb_torsion.get_atom_ids()
                    ]
                except KeyError:
                    # The atoms of the building block torsion do not
                    # all have corresponding atoms in the constructed
                    # molecule (e.g. atom was deleted in construction)
                    # so skip to next building block torsion.
                    continue
                matched_atom_ids = [
                    atom.get_id()
                    for atom in matched_atoms
                ]
                if set(matched_atom_ids[1:3]) == set(atom_ids[1:3]):
                    # Set the constructed molecule torsion to match the
                    # building block torsion.
                    torsions[i] = Torsion(*matched_atoms)
                    break

        yield tuple(torsions)
