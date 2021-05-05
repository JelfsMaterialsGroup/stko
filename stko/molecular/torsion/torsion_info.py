"""
TorsionInfo
===========

Class defining torsion angle info.

"""

import logging

logger = logging.getLogger(__name__)


class TorsionInfo:
    """
    Holds additional info about ConstructedMoleculeTorsioned torsions.

    Attributes
    ----------
    torsion:
        The torsion about which information is held.

    building_block_torsion:
        The building block torsion from which this torsion originates.
        Can be ``None``, if the atoms that make up the torsion did not
        come from a single building block.

    building_block:
        The building block from which this torsion originates.
        Can be ``None``, if the atoms that make up the torsion did not
        come from a single building block.

    building_block_id:
        A unique id for each :class:`.Molecule` placed during
        the construction of the :class:`.ConstructedMolecule`. As a
        single :class:`.Molecule` can be placed multiple times
        during construction, the `building_block_id` allows
        the user to distinguish between each placement. Can be
        Can be ``None``, if the atoms that make up the torsion did not
        come from a single building block.

    """

    def __init__(
        self,
        torsion,
        building_block,
        building_block_id,
        building_block_torsion,
    ):
        self._torsion = torsion
        self._building_block = building_block
        self._building_block_id = building_block_id
        self._building_block_torsion = building_block_torsion

    def get_torsion(self):
        """
        Torsion of atoms in constructed molecule.

        """

        return self._torsion

    def get_building_block(self):
        return self._building_block

    def get_building_block_torsion(self):
        """
        Torsion of atoms in building block, not constructed mol.

        """

        return self._building_block_torsion

    def get_building_block_id(self):
        return self._building_block_id
