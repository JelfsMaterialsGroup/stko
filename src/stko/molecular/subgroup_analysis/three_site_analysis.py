import logging

import stk
import numpy as np

from stko.functional_groups import ThreeSiteFG
from stko.utilities.exceptions import NotDitopicThreeSiteError
from stko.utilities.utilities import (
    get_atom_distance,
    vector_angle,
)

logger = logging.getLogger(__name__)


class DitopicThreeSiteAnalyser:
    """
    Analyses geometry of functional groups in ditopic molecules.

    """

    def _check_functional_groups(self, molecule: stk.BuildingBlock) -> None:
        """
        Check if the molecule has two ditopic functional groups.

        Raises:

            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """

        fg_counts = 0
        for fg in molecule.get_functional_groups():
            if isinstance(fg, ThreeSiteFG):  # type: ignore[unreachable]
                fg_counts += 1  # type: ignore[unreachable]

        if fg_counts != 2:
            raise NotDitopicThreeSiteError(
                f"{molecule} does not have 2 ThreeSiteFG functional groups."
            )

    def get_binder_distance(self, molecule: stk.BuildingBlock) -> float:
        """
        Get the distance between binder atoms in Angstrom.

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
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Get the position of centroids of atoms adjacent to binder atoms.

        Parameters:

            molecule:
                Molecule to analyse.

        Raises:

            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """

        self._check_functional_groups(molecule)

        adj_centroids = []
        for fg in molecule.get_functional_groups():
            adj_atom_ids = (fg.get_neigh1().get_id(), fg.get_neigh2().get_id())
            adj_centroids.append(molecule.get_centroid(atom_ids=adj_atom_ids))
        return tuple(adj_centroids)

    def get_binder_centroid_angle(
        self,
        molecule: stk.BuildingBlock,
    ) -> float:
        """
        Get the angle between binders and molecule centroid.

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
        angle = np.degrees(
            vector_angle(vector1=fg_vectors[0], vector2=fg_vectors[1])
        )
        return angle

    def get_binder_binder_angle(self, molecule: stk.BuildingBlock) -> float:
        """
        Get the angle between binders-adjacent vectors.

        Parameters:

            molecule:
                Molecule to analyse.

        Raises:

            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """

        self._check_functional_groups(molecule)

        raise NotImplementedError()

    def get_binder_angles(
        self,
        molecule: stk.BuildingBlock,
    ) -> tuple[float, float]:
        """
        Get the binder-adjacent-binder angles.

        Represents the two `reaction angle`s of the ligand.

        Parameters:

            molecule:
                Molecule to analyse.

        Raises:

            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """

        self._check_functional_groups(molecule)

        raise NotImplementedError()

    def get_binder_adjacent_torsion(
        self, molecule: stk.BuildingBlock
    ) -> float:
        """
        Get the torsion between binders and their adjacent centroids.

        Parameters:

            molecule:
                Molecule to analyse.

        Raises:

            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """

        self._check_functional_groups(molecule)

        raise NotImplementedError()
