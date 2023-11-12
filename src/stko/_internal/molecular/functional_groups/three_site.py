import logging
from collections import abc
from typing import Self

import stk

logger = logging.getLogger(__name__)


class ThreeSiteFG:
    """
    Represents FG sites like N atom in pyridine functional group.

    WARNING: This code is only present in the latest versions of stko
    that require Python 3.11!

    The structure of the functional group is given by the pseudo-SMILES
    ``[neighbour][binder][neighbour]``.

    Contains :class:`stk.GenericFunctionalGroup`.

    """

    def __init__(
        self,
        neigh1: stk.Atom,
        binder: stk.Atom,
        neigh2: stk.Atom,
        bonders: tuple[stk.Atom, ...],
        deleters: tuple[stk.Atom, ...],
    ) -> None:
        """
        Parameters:

            neigh1:
                The first neighbour atom.

            binder:
                The central atom that forms the bond.

            neigh2:
                The second neighbour atom.

            bonders:
                The bonder atoms, this should be the same ID as `binder`.

            deleters:
                The deleter atoms, there should be none.

        """

        self._neigh1 = neigh1
        self._binder = binder
        self._neigh2 = neigh2
        atoms = (neigh1, binder, neigh2)
        self._functional_group = stk.GenericFunctionalGroup(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
        )

    def get_neigh1(self) -> stk.Atom:
        return self._neigh1

    def get_neigh2(self) -> stk.Atom:
        return self._neigh2

    def get_binder(self) -> stk.Atom:
        return self._binder

    def clone(self) -> Self:
        clone = self.__class__.__new__(self.__class__)
        clone._neigh1 = self._neigh1
        clone._binder = self._binder
        clone._neigh2 = self._neigh2
        clone._functional_group = self._functional_group.clone()
        return clone

    def get_bonders(self) -> abc.Iterator[stk.Atom]:
        yield from self._functional_group._bonders

    def get_num_bonders(self) -> int:
        return len(self._functional_group._bonders)

    def get_bonder_ids(self) -> abc.Iterator[int]:
        yield from (a.get_id() for a in self._functional_group._bonders)

    def get_deleters(self) -> abc.Iterator[stk.Atom]:
        yield from self._functional_group._deleters

    def get_deleter_ids(self) -> abc.Iterator[int]:
        yield from (a.get_id() for a in self._functional_group._deleters)

    def get_atoms(self) -> abc.Iterator[stk.Atom]:
        yield from self._functional_group._atoms

    def get_atom_ids(self) -> abc.Iterator[int]:
        yield from (a.get_id() for a in self._functional_group._atoms)

    def get_placer_ids(self) -> abc.Iterator[int]:
        yield from (a.get_id() for a in self._functional_group._placers)

    def get_core_atom_ids(self) -> abc.Iterator[int]:
        yield from (a.get_id() for a in self._functional_group._core_atoms)

    def with_atoms(self, atom_map: dict[int, stk.Atom]) -> Self:
        clone = self.__class__.__new__(self.__class__)
        clone._functional_group = stk.GenericFunctionalGroup(
            atoms=tuple(
                atom_map.get(a.get_id(), a)
                for a in self._functional_group._atoms
            ),
            bonders=tuple(
                atom_map.get(a.get_id(), a)
                for a in self._functional_group._bonders
            ),
            deleters=tuple(
                atom_map.get(a.get_id(), a)
                for a in self._functional_group._deleters
            ),
            placers=tuple(
                atom_map.get(a.get_id(), a)
                for a in self._functional_group._placers
            ),
        )
        clone._neigh1 = atom_map.get(self._neigh1.get_id(), self._neigh1)
        clone._binder = atom_map.get(self._binder.get_id(), self._binder)
        clone._neigh2 = atom_map.get(self._neigh2.get_id(), self._neigh2)
        return clone

    def with_ids(self, id_map: dict[int, int]) -> Self:
        clone = self.__class__.__new__(self.__class__)
        clone._functional_group = self._functional_group.with_ids(id_map)
        clone._neigh1 = self._neigh1.with_id(
            id_map.get(self._neigh1.get_id(), self._neigh1.get_id())
        )
        clone._binder = self._binder.with_id(
            id_map.get(self._binder.get_id(), self._binder.get_id())
        )
        clone._neigh2 = self._neigh2.with_id(
            id_map.get(self._neigh2.get_id(), self._neigh2.get_id())
        )
        return clone

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"{self._neigh1}, {self._binder}, {self._neigh2}, "
            f"bonders={self._functional_group._bonders})"
        )


class ThreeSiteFactory(stk.FunctionalGroupFactory):
    """
    Find ThreeSite functional groups in molecules.

    WARNING: This code is only present in the latest versions of stko
    that require Python 3.11!


    """

    def __init__(
        self,
        smarts: str,
        bonders: tuple[int, ...] = (1,),
        deleters: tuple[int, ...] = (),
    ) -> None:
        """
        Parameters:

            smarts:
                SMARTS string to use to find functional group. Of form
                ``[neighbour][binder][neighbour]``.

            bonders:
                The bonder atoms, this should be the same ID as `binder`.

            deleters:
                The deleter atoms, there should be none.

        """

        self._smarts = smarts
        self._bonders = bonders
        self._deleters = deleters

    def get_functional_groups(
        self,
        molecule: stk.Molecule,
    ) -> abc.Iterable[ThreeSiteFG]:
        generic_functional_groups = stk.SmartsFunctionalGroupFactory(
            smarts=self._smarts,
            bonders=self._bonders,
            deleters=self._deleters,
        ).get_functional_groups(molecule)
        for fg in generic_functional_groups:
            atom_ids = (i.get_id() for i in fg.get_atoms())
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield ThreeSiteFG(
                neigh1=atoms[0],
                binder=atoms[1],
                neigh2=atoms[2],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )


class CNCFactory(ThreeSiteFactory):
    """
    A subclass of :class:`.ThreeSiteFactory`.

    WARNING: This code is only present in the latest versions of stko
    that require Python 3.11!

    SMARTs string for [carbon][nitrogen][carbon]: "[#6]~[#7X2]~[#6]"

    """

    def __init__(
        self,
        bonders: tuple[int, ...] = (1,),
        deleters: tuple[int, ...] = (),
    ):
        self._smarts = "[#6]~[#7X2]~[#6]"
        self._bonders = bonders
        self._deleters = deleters


class CNNFactory(ThreeSiteFactory):
    """
    A subclass of :class:`.ThreeSiteFactory`.

    WARNING: This code is only present in the latest versions of stko
    that require Python 3.11!

    SMARTs string for [nitrogen][nitrogen][carbon]: "[#7]~[#7X2]~[#6]"

    """

    def __init__(
        self,
        bonders: tuple[int, ...] = (1,),
        deleters: tuple[int, ...] = (),
    ):
        self._smarts = "[#7]~[#7X2]~[#6]"
        self._bonders = bonders
        self._deleters = deleters


class NNNFactory(ThreeSiteFactory):
    """
    A subclass of :class:`.ThreeSiteFactory`.

    WARNING: This code is only present in the latest versions of stko
    that require Python 3.11!

    SMARTs string for [nitrogen][nitrogen][nitrogen]: "[#7]~[#7X2]~[#7]"

    """

    def __init__(
        self,
        bonders: tuple[int, ...] = (1,),
        deleters: tuple[int, ...] = (),
    ):
        self._smarts = "[#7]~[#7X2]~[#7]"
        self._bonders = bonders
        self._deleters = deleters
