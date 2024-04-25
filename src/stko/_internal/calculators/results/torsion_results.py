from collections import abc, defaultdict

import stk

from stko._internal.molecular.torsion.torsion import Torsion
from stko._internal.molecular.torsion.torsion_info import TorsionInfo
from stko._internal.utilities.utilities import calculate_dihedral


class TorsionResults:
    """Results class containing molecule torsions."""

    def __init__(self, generator: abc.Generator, mol: stk.Molecule) -> None:
        self._torsions = next(generator)
        self._mol = mol

    def get_torsions(self) -> abc.Iterable[Torsion]:
        return self._torsions

    def get_molecule(self) -> stk.Molecule:
        return self._mol

    def get_torsion_angles(self) -> abc.Iterable[tuple[Torsion, float]]:
        for torsion in self._torsions:
            yield (
                torsion,
                calculate_dihedral(
                    pt1=next(
                        self._mol.get_atomic_positions(
                            torsion.get_atom_ids()[0]
                        )
                    ),
                    pt2=next(
                        self._mol.get_atomic_positions(
                            torsion.get_atom_ids()[1]
                        )
                    ),
                    pt3=next(
                        self._mol.get_atomic_positions(
                            torsion.get_atom_ids()[2]
                        )
                    ),
                    pt4=next(
                        self._mol.get_atomic_positions(
                            torsion.get_atom_ids()[3]
                        )
                    ),
                ),
            )


class ConstructedMoleculeTorsionResults(TorsionResults):
    """Results class containing molecule torsions."""

    def __init__(
        self,
        generator: abc.Generator,
        mol: stk.ConstructedMolecule,
    ) -> None:
        self._torsions = next(generator)
        self._mol: stk.ConstructedMolecule = mol

    def get_torsion_infos_by_building_block(
        self,
    ) -> dict[int | None, list[TorsionInfo]]:
        """Returns dictionary of torsions by building block."""
        torsion_infos_by_building_block = defaultdict(list)
        for torsion_info in self.get_torsion_infos():
            if torsion_info.get_building_block_id() is not None:
                torsion_infos_by_building_block[
                    torsion_info.get_building_block_id()
                ].append(torsion_info)
        return torsion_infos_by_building_block

    def get_torsion_infos(self) -> abc.Iterable[TorsionInfo]:
        for torsion in self._torsions:
            atom_infos = list(
                self._mol.get_atom_infos(
                    atom_ids=(i for i in torsion.get_atom_ids())
                )
            )
            # Get atom info and check they are all the same.
            building_block_ids = {
                i.get_building_block_id() for i in atom_infos
            }
            if len(building_block_ids) > 1:
                same_building_block = False
            else:
                same_building_block = True

            if same_building_block:
                building_block_id = next(iter(building_block_ids))

                building_block = next(
                    i.get_building_block() for i in atom_infos
                )
                bb_atoms = tuple(
                    i.get_building_block_atom() for i in atom_infos
                )
                building_block_torsion = Torsion(
                    bb_atoms[0], bb_atoms[1], bb_atoms[2], bb_atoms[3]
                )
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
