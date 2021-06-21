"""
Torsion Calculators
===================

#. :class:`.TorsionCalculator`
#. :class:`.ConstructedMoleculeTorsionCalculator`

Methods to extract torsions from a molecule or constructed molecule.

"""

from collections import defaultdict
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

        # Create a molecule whose torsions we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the calculator.
        tc = stko.TorsionCalculator()

        # Extract the torsions.
        tc_results = tc.get_results(mol1)
        for t, ang in tc_results.get_torsion_angles():
            print(t, ang, t.get_atom_ids())

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

    """
    
    def calculate(self, mol):
        """extract torsions with rdkit, then match to building blocks
        """
        def get_atom_maps():
            """
            map from building block atom ids to constructed molecule atoms for
            a specified building block id
            """
            atom_maps = defaultdict(dict)
            for atom_info in mol.get_atom_infos():
                current_atom_map = atom_maps[atom_info.get_building_block_id()]
                bb_atom_id = atom_info.get_building_block_atom().get_id()
                current_atom_map[bb_atom_id] = atom_info.get_atom()
            return atom_maps
        
        torsions = list(next(super().calculate(mol)))
        atom_maps = get_atom_maps()
        
        # loop over torsions, updating each to match a building block
        # torsion if possible
        for i, torsion in enumerate(torsions):
            # atom ids, atom infos, and building block ids of current
            # torsion in constructed molecule
            atom_ids = list(torsion.get_atom_ids())
            atom_infos = list(mol.get_atom_infos(atom_ids))
            build_block_ids = [atom_info.get_building_block_id()
                               for atom_info in atom_infos]
            
            # check if two central atoms of torsion are from the same building
            # block. if not, leave this torsion alone
            if build_block_ids[1] is None:
                continue
            if (atom_infos[1].get_building_block_id()
                    != atom_infos[2].get_building_block_id()):
                continue
            
            central_atom_ids = set(atom_ids[1:3])
            build_block_torsions = TorsionCalculator().get_results(
                    atom_infos[1].get_building_block()).get_torsions()
            
            # look for a torsion in the building block that has the same
            # central atoms
            for bb_torsion in build_block_torsions:
                bb_central_atom_ids = set(list(bb_torsion.get_atom_ids())[1:3])
                if bb_central_atom_ids == central_atom_ids:
                    # set the constructed molecule torsion to match the
                    # building block torsion
                    atoms = [atom_maps[build_block_ids[1]][atom_id]
                             for atom_id in bb_torsion.get_atom_ids()]
                    torsions[i] = Torsion(*atoms)
                    break
            
        yield tuple(torsions)

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
