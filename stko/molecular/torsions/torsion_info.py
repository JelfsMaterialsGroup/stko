from dataclasses import dataclass
from stko.molecular.torsions.torsion import Torsion
import stk
@dataclass
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
    torsion: Torsion
    building_block_torsion: Torsion
    building_block: stk.BuildingBlock
    building_block_id: int
    