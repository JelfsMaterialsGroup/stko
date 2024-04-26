import logging
from collections import abc

import numpy as np
import stk

from stko._internal.calculators.planarity_calculators import (
    PlanarityCalculator,
)
from stko._internal.utilities.utilities import vector_angle

logger = logging.getLogger(__name__)


class Subgroup:
    """Defines subgroup to search for and its measure of interest.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def _find_subgroup(
        self,
        molecule: stk.Molecule,
    ) -> abc.Iterable[stk.FunctionalGroup]:
        raise NotImplementedError

    def _calculate_measure(
        self,
        molecule: stk.Molecule,
        atom_group: stk.FunctionalGroup,
    ) -> float:
        raise NotImplementedError

    def measure(self, molecule: stk.Molecule) -> list[float]:
        """Measure the geometrical property of interest for the subgroup.

        Defined in `_calculate_measure` and `_find_subgroup` private
        methods.

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            List of measures.

        """
        atom_groups = self._find_subgroup(molecule)
        return [self._calculate_measure(molecule, ag) for ag in atom_groups]


class C6Planarity(Subgroup):
    """Sub group measure of C6 ring planarity.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def _find_subgroup(
        self,
        molecule: stk.Molecule,
    ) -> abc.Iterable[stk.FunctionalGroup]:
        smarts = "[#6X3]1@[#6X3]@[#6X3]@[#6X3]@[#6X3]@[#6X3]1"
        with_fgs = stk.BuildingBlock.init_from_molecule(
            molecule=molecule,
            functional_groups=(
                stk.SmartsFunctionalGroupFactory(
                    smarts=smarts,
                    bonders=(),
                    deleters=(),
                ),
            ),
        )
        yield from with_fgs.get_functional_groups()

    def _calculate_measure(
        self,
        molecule: stk.Molecule,
        atom_group: stk.FunctionalGroup,
    ) -> float:
        atom_ids = [i.get_id() for i in atom_group.get_atoms()]
        # Create the calculator.
        pc = PlanarityCalculator()
        # Extract the measures.
        pc_results = pc.get_results(
            molecule,
            plane_atom_ids=atom_ids,
            deviation_atom_ids=atom_ids,
        )
        return pc_results.get_plane_deviation()


class C5N1Planarity(C6Planarity):
    """Sub group measure of C5N1 ring planarity.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def _find_subgroup(
        self,
        molecule: stk.Molecule,
    ) -> abc.Iterable[stk.FunctionalGroup]:
        smarts = "[#6X3]1@[#7X2]@[#6X3]@[#6X3]@[#6X3]@[#6X3]1"
        with_fgs = stk.BuildingBlock.init_from_molecule(
            molecule=molecule,
            functional_groups=(
                stk.SmartsFunctionalGroupFactory(
                    smarts=smarts,
                    bonders=(),
                    deleters=(),
                ),
            ),
        )
        yield from with_fgs.get_functional_groups()


class X5Planarity(C6Planarity):
    """Sub group measure of five-membered ring planarity.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def _find_subgroup(
        self,
        molecule: stk.Molecule,
    ) -> abc.Iterable[stk.FunctionalGroup]:
        smarts = "[*]1@[*]@[*]@[*]@[*]1"
        with_fgs = stk.BuildingBlock.init_from_molecule(
            molecule=molecule,
            functional_groups=(
                stk.SmartsFunctionalGroupFactory(
                    smarts=smarts,
                    bonders=(),
                    deleters=(),
                ),
            ),
        )
        yield from with_fgs.get_functional_groups()


class AlkyneAngle(Subgroup):
    """Sub group measure of avg. CCC angle in C-C#C-C group.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def _find_subgroup(
        self,
        molecule: stk.Molecule,
    ) -> abc.Iterable[stk.FunctionalGroup]:
        smarts = "[#6][#6]#[#6][#6]"
        with_fgs = stk.BuildingBlock.init_from_molecule(
            molecule=molecule,
            functional_groups=(
                stk.SmartsFunctionalGroupFactory(
                    smarts=smarts,
                    bonders=(),
                    deleters=(),
                ),
            ),
        )
        yield from with_fgs.get_functional_groups()

    def _calculate_measure(
        self,
        molecule: stk.Molecule,
        atom_group: stk.FunctionalGroup,
    ) -> float:
        atom_ids = [i.get_id() for i in atom_group.get_atoms()]
        atom_positions = [
            next(molecule.get_atomic_positions(atom_ids=i)) for i in atom_ids
        ]
        vector_sets = [
            (
                atom_positions[0] - atom_positions[1],
                atom_positions[2] - atom_positions[1],
            ),
            (
                atom_positions[1] - atom_positions[2],
                atom_positions[3] - atom_positions[2],
            ),
        ]

        return np.average(
            [
                np.degrees(vector_angle(vects[0], vects[1]))
                for vects in vector_sets
            ]
        )


class SubgroupAnalyser:
    """Defines the analyser of subgroup instances.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def calculate(self, molecule: stk.Molecule) -> dict[str, list[float]]:
        """Calculate the all possible fg values in molecule."""
        fg_definitions: dict[str, Subgroup] = {
            "c6_planarity": C6Planarity(),
            "c5n1_planarity": C5N1Planarity(),
            "x5_planarity": X5Planarity(),
            "c#c_angle": AlkyneAngle(),
        }

        fg_results = {}
        for pot_fg in fg_definitions:
            fg_cls = fg_definitions[pot_fg]
            fg_results[pot_fg] = fg_cls.measure(molecule)

        return fg_results
