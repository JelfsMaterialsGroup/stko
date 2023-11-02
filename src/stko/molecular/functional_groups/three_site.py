import logging

import stk

logger = logging.getLogger(__name__)


class ThreeSiteFG(stk.GenericFunctionalGroup):
    """
    Represents FG sites like N atom in pyridine functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[neighbour][binder][neighbour]``.

    """

    def __init__(
        self,
        neigh1: stk.Atom,
        binder: stk.Atom,
        neigh2: stk.Atom,
        bonders: tuple[stk.Atom],
        deleters: tuple[stk.Atom],
    ):
        """
        Initialize a :class:`.ThreeSiteFG` instance.

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
        super().__init__(atoms, bonders, deleters)

    def get_neigh1(self):
        return self._neigh1

    def get_neigh2(self):
        return self._neigh2

    def get_binder(self):
        return self._binder

    def clone(self):
        clone = super().clone()
        clone._neigh1 = self._neigh1
        clone._binder = self._binder
        clone._neigh2 = self._neigh2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._neigh1 = atom_map.get(
            self._neigh1.get_id(),
            self._neigh1,
        )
        clone._nitrogen = atom_map.get(
            self._binder.get_id(),
            self._binder,
        )
        clone._neigh2 = atom_map.get(
            self._neigh2.get_id(),
            self._neigh2,
        )
        return clone

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            f"{self._neigh1}, {self._binder}, {self._neigh2}, "
            f"bonders={self._bonders})"
        )


class ThreeSiteFactory(stk.FunctionalGroupFactory):
    """
    A subclass of :class:`stk.FunctionalGroupFactory`.

    """

    def __init__(
        self,
        smarts: str,
        bonders: tuple[int] = (1,),
        deleters: tuple[int] = (),
    ) -> None:
        """
        Intiailize :class:`ThreeSiteFactory`.

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

    def get_functional_groups(self, molecule):
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

    SMARTs string for [carbon][nitrogen][carbon]: "[#6]~[#7X2]~[#6]"

    """

    def __init__(self, bonders: tuple[int] = (1,), deleters: tuple[int] = ()):
        """
        Intiailize :class:`CNCFactory`.

        """
        self._smarts = "[#6]~[#7X2]~[#6]"
        self._bonders = bonders
        self._deleters = deleters


class CNNFactory(ThreeSiteFactory):
    """
    A subclass of :class:`.ThreeSiteFactory`.

    SMARTs string for [nitrogen][nitrogen][carbon]: "[#7]~[#7X2]~[#6]"

    """

    def __init__(self, bonders: tuple[int] = (1,), deleters: tuple[int] = ()):
        """
        Intiailize :class:`CNNFactory`.

        """
        self._smarts = "[#7]~[#7X2]~[#6]"
        self._bonders = bonders
        self._deleters = deleters


class NNNFactory(ThreeSiteFactory):
    """
    A subclass of :class:`.ThreeSiteFactory`.

    SMARTs string for [nitrogen][nitrogen][nitrogen]: "[#7]~[#7X2]~[#7]"

    """

    def __init__(self, bonders: tuple[int] = (1,), deleters: tuple[int] = ()):
        """
        Intiailize :class:`NNNFactory`.

        """
        self._smarts = "[#7]~[#7X2]~[#7]"
        self._bonders = bonders
        self._deleters = deleters
