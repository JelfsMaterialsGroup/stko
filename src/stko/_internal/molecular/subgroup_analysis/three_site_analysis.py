import logging

import numpy as np
import stk

from stko._internal.molecular.functional_groups.three_site import ThreeSiteFG
from stko._internal.utilities.exceptions import NotDitopicThreeSiteError
from stko._internal.utilities.utilities import (
    calculate_dihedral,
    get_atom_distance,
    vector_angle,
)

logger = logging.getLogger(__name__)


class DitopicThreeSiteAnalyser:
    """Analyses geometry of functional groups in ditopic molecules.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def _check_functional_groups(self, molecule: stk.BuildingBlock) -> None:
        """Check if the molecule has two ditopic functional groups.

        Raises:
            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """
        fg_counts = 0
        for fg in molecule.get_functional_groups():
            if isinstance(fg, ThreeSiteFG):  # type: ignore[unreachable]
                fg_counts += 1  # type: ignore[unreachable]

        if fg_counts != 2:  # noqa: PLR2004
            msg = f"{molecule} does not have 2 ThreeSiteFG functional groups."
            raise NotDitopicThreeSiteError(msg)

    def get_binder_distance(self, molecule: stk.BuildingBlock) -> float:
        """Get the distance between binder atoms in Angstrom.

        Parameters:
            molecule:
                Molecule to analyse.

        Raises:
            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """
        self._check_functional_groups(molecule)

        binder_ids = [
            fg.get_binder().get_id()  # type: ignore[attr-defined]
            for fg in molecule.get_functional_groups()
        ]
        return get_atom_distance(
            position_matrix=molecule.get_position_matrix(),
            atom1_id=binder_ids[0],
            atom2_id=binder_ids[1],
        )

    def get_adjacent_centroids(
        self,
        molecule: stk.BuildingBlock,
    ) -> tuple[np.ndarray, ...]:
        """Get the position of centroids of atoms adjacent to binder atoms.

        Parameters:
            molecule:
                Molecule to analyse.

        Raises:
            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """
        self._check_functional_groups(molecule)

        adj_centroids = [
            self._get_adj_centroid(molecule, fg)
            for fg in molecule.get_functional_groups()
        ]
        return tuple(adj_centroids)

    def get_binder_centroid_angle(
        self,
        molecule: stk.BuildingBlock,
    ) -> float:
        """Get the angle between binders and molecule centroid.

        Parameters:
            molecule:
                Molecule to analyse.

        Raises:
            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """
        self._check_functional_groups(molecule)

        position_matrix = molecule.get_position_matrix()
        binder_ids = [
            fg.get_binder().get_id()  # type: ignore[attr-defined]
            for fg in molecule.get_functional_groups()
        ]
        binder_positions = [position_matrix[id_] for id_ in binder_ids]

        # Get building block centroid.
        centroid_position = molecule.get_centroid()

        # Get vectors.
        fg_vectors = [i - centroid_position for i in binder_positions]

        # Calculate the angle between the two vectors.
        return np.degrees(
            vector_angle(vector1=fg_vectors[0], vector2=fg_vectors[1])
        )

    def get_binder_binder_angle(self, molecule: stk.BuildingBlock) -> float:
        """Get the angle between binders-adjacent vectors.

        Parameters:
            molecule:
                Molecule to analyse.

        Raises:
            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """
        self._check_functional_groups(molecule)

        vectors = [
            self._get_binder_adjacent_vector(molecule, fg, normalised=True)
            for fg in molecule.get_functional_groups()
        ]
        return np.degrees(vector_angle(vector1=vectors[0], vector2=vectors[1]))

    def get_binder_angles(
        self,
        molecule: stk.BuildingBlock,
    ) -> tuple[float, float]:
        """Get the binder-adjacent-binder angles.

        Represents the two reaction angles of the ligand, or the internal
        angles from DOI: 10.1039/D3SC03991A.

        Caution with the direction of the vectors!

        Parameters:
            molecule:
                Molecule to analyse.

        Raises:
            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """
        self._check_functional_groups(molecule)

        binder_centroids = [
            self._get_binder_centroid(molecule, fg)
            for fg in molecule.get_functional_groups()
        ]
        binder_vector = binder_centroids[1] - binder_centroids[0]

        vectors = [
            self._get_binder_adjacent_vector(molecule, fg, normalised=False)
            for fg in molecule.get_functional_groups()
        ]
        return (
            180 - np.degrees(vector_angle(vectors[0], -binder_vector)),
            180 - np.degrees(vector_angle(vectors[1], binder_vector)),
        )

    def get_halfbite_angles(
        self,
        molecule: stk.BuildingBlock,
    ) -> tuple[float, float]:
        """Get the half-bite angles defined by the binders.

        The bite angle is a common measure used in metal-organic cages
        that represents the angle between the two reaction angles of
        the ligand as if the molecule was symmetric and the torsion
        between binders is 0 degrees. Caution using this for many
        molecules that are flexible!

        Here the bite angle is defined ``visually`` as in:

            https://doi.org/10.1016/j.ccr.2018.06.010

        Technically, the measure of interest is the sum of the two
        floats output by this method, which should be two times either
        one. (They should be similar!)

        Parameters:
            molecule:
                Molecule to analyse.

        Raises:
            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """
        binder_angles = self.get_binder_angles(molecule)
        return (binder_angles[0] - 90, binder_angles[1] - 90)

    def get_binder_adjacent_torsion(
        self, molecule: stk.BuildingBlock
    ) -> float:
        """Get the torsion (-180, 180) between binders and their adjacents.

        Parameters:
            molecule:
                Molecule to analyse.

        Raises:
            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """
        self._check_functional_groups(molecule)

        binder_centroids = []
        adj_centroids = []
        for fg in molecule.get_functional_groups():
            binder_centroids.append(self._get_binder_centroid(molecule, fg))
            adj_centroids.append(self._get_adj_centroid(molecule, fg))

        return calculate_dihedral(
            pt1=binder_centroids[0],
            pt2=adj_centroids[0],
            pt3=adj_centroids[1],
            pt4=binder_centroids[1],
        )

    def _get_binder_centroid(
        self,
        molecule: stk.Molecule,
        fg: stk.FunctionalGroup,
    ) -> np.ndarray:
        return molecule.get_centroid(
            fg.get_binder().get_id()  # type: ignore[attr-defined]
        )

    def _get_adj_centroid(
        self,
        molecule: stk.Molecule,
        fg: stk.FunctionalGroup,
    ) -> np.ndarray:
        adj_atom_ids = (
            fg.get_neigh1().get_id(),  # type: ignore[attr-defined]
            fg.get_neigh2().get_id(),  # type: ignore[attr-defined]
        )
        return molecule.get_centroid(atom_ids=adj_atom_ids)

    def _get_binder_adjacent_vector(
        self,
        molecule: stk.BuildingBlock,
        fg: stk.FunctionalGroup,
        normalised: bool = True,
    ) -> np.ndarray:
        binder_centroid = self._get_binder_centroid(molecule, fg)
        adj_centroid = self._get_adj_centroid(molecule, fg)
        if normalised:
            return (binder_centroid - adj_centroid) / np.linalg.norm(
                binder_centroid - adj_centroid
            )

        return binder_centroid - adj_centroid
