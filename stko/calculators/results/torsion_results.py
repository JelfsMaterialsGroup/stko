"""
Torsion Results
===============

#. :class:`.TorsionResults`
#. :class:`.ConstructedMoleculeTorsionResults`

Results classes for extracting molecular torsions.

"""

from collections import defaultdict
from .results import Results
from ...molecular.torsion import TorsionInfo, Torsion
from ...utilities import calculate_dihedral


class TorsionResults(Results):
    """
    Results class containing molecule torsions.

    """

    def __init__(self, generator, mol):
        self._torsions = next(generator)
        self._mol = mol

    def get_torsions(self):
        return self._torsions

    def get_molecule(self):
        return self._mol

    def get_torsion_angles(self):
        for torsion in self._torsions:
            print('a', torsion)
            yield (
                torsion, calculate_dihedral(
                    pt1=tuple(
                        self._mol.get_atomic_positions(
                            torsion.get_atom_ids()[0]
                        )
                    )[0],
                    pt2=tuple(
                        self._mol.get_atomic_positions(
                            torsion.get_atom_ids()[1]
                        )
                    )[0],
                    pt3=tuple(
                        self._mol.get_atomic_positions(
                            torsion.get_atom_ids()[2]
                        )
                    )[0],
                    pt4=tuple(
                        self._mol.get_atomic_positions(
                            torsion.get_atom_ids()[3]
                        )
                    )[0],
                )
            )


class ConstructedMoleculeTorsionResults(TorsionResults):
    """
    Results class containing molecule torsions.

    """

    def __init__(self, generator, mol):
        self._torsions = next(generator)
        self._mol = mol

    def get_torsion_infos_by_building_block(self):
        """
        Returns dictionary of torsions by building block.

        """

        torsion_infos_by_building_block = defaultdict(list)
        for torsion_info in self.get_torsion_infos():
            if torsion_info.get_building_block_id() is not None:
                torsion_infos_by_building_block[
                    torsion_info.get_building_block_id()
                ].append(torsion_info)
        return torsion_infos_by_building_block

    def get_torsion_infos(self):
        for torsion in self._torsions:
            atom_infos = list(
                self._mol.get_atom_infos(
                    atom_ids=(i for i in torsion.get_atom_ids())
                )
            )
            # Get atom info and check they are all the same.
            building_block_ids = set((
                i.get_building_block_id() for i in atom_infos
            ))
            if len(building_block_ids) > 1:
                same_building_block = False
            else:
                same_building_block = True

            if same_building_block:
                building_block_id = next(iter(building_block_ids))

                building_block = tuple((
                   i.get_building_block() for i in atom_infos
                ))[0]
                bb_atoms = tuple(
                    i.get_building_block_atom()
                    for i in atom_infos
                )
                building_block_torsion = Torsion(*bb_atoms)
                yield TorsionInfo(
                    torsion=torsion,
                    building_block=building_block,
                    building_block_id=building_block_id,
                    building_block_torsion=building_block_torsion,
                )
            else:
                yield TorsionInfo(
                    torsion=torsion,
                    building_block=None,
                    building_block_id=None,
                    building_block_torsion=None,
                )
