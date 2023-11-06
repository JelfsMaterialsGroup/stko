import itertools as it
import logging

import numpy as np
import stk
from scipy.spatial.distance import cdist

from stko.utilities.utilities import get_atom_distance, vector_angle

logger = logging.getLogger(__name__)


class PoreAnalyser:
    """
    Analyses geometry of building blocks in constructed molecules.

    """

    def _get_metal_atom_ids(
        self,
        molecule: stk.Molecule,
        metal_atom_nos: tuple[int],
    ) -> list[int]:
        return [
            i.get_id()
            for i in molecule.get_atoms()
            if i.get_atomic_number() in metal_atom_nos
        ]

    def get_metal_distances(
        self,
        molecule: stk.Molecule,
        metal_atom_nos: tuple[int],
    ) -> dict[tuple[int, int], float]:
        """
        Get all metal atom pair distances.

        Parameters:

            molecule:
                The molecule to analyse.

            metal_atom_nos:
                The atomic numbers to delete. Can be a tuple of one or
                any element on periodic table.

        Returns:

            The distances and associated metal atom ids.

        """

        metal_atom_ids = self._get_metal_atom_ids(molecule, metal_atom_nos)
        position_matrix = molecule.get_position_matrix()

        distances = {}
        for a1id, a2id in it.combinations(metal_atom_ids, 2):
            distances[(a1id, a2id)] = get_atom_distance(
                position_matrix=position_matrix,
                atom1_id=a1id,
                atom2_id=a2id,
            )

        return distances

    def get_metal_centroid_metal_angle(
        self,
        molecule: stk.Molecule,
        metal_atom_nos: tuple[int],
    ) -> dict[tuple[int, int], float]:
        """
        Get all metal-centroid-metal angles.

        Parameters:

            molecule:
                The molecule to analyse.

            metal_atom_nos:
                The atomic numbers to delete. Can be a tuple of one or
                any element on periodic table.

        Returns:

            The angles in degrees and associated metal atom ids.

        """

        metal_atom_ids = self._get_metal_atom_ids(molecule, metal_atom_nos)
        position_matrix = molecule.get_position_matrix()
        centroid = molecule.get_centroid()

        angles = {}
        for a1id, a2id in it.combinations(metal_atom_ids, 2):
            vector1 = position_matrix[a1id] - centroid
            vector2 = position_matrix[a2id] - centroid
            angles[(a1id, a2id)] = np.degrees(vector_angle(vector1, vector2))

        return angles

    def get_min_centroid_distance(
        self,
        molecule: stk.Molecule,
    ) -> float:
        """
        Get the minimum distance between the molecule and centroid.

        This is nearly equivalent to a pore radius.

        Parameters:

            molecule:
                The molecule to analyse.

        Returns:

            The minimum centroid to atom distance.

        """

        pair_dists = cdist(
            molecule.get_position_matrix(),
            molecule.get_centroid().reshape(1, 3),
        )

        return np.min(pair_dists.flatten())

    def get_avg_centroid_distance(
        self,
        molecule: stk.Molecule,
    ) -> tuple[float, float]:
        """
        Get the average distance between the molecule and centroid.

        Parameters:

            molecule:
                The molecule to analyse.

        Returns:

            The average and std. deviation of centroid to atom distances.

        """

        pair_dists = cdist(
            molecule.get_position_matrix(),
            molecule.get_centroid().reshape(1, 3),
        )
        flattened = pair_dists.flatten()

        return (np.mean(flattened), np.std(flattened))
