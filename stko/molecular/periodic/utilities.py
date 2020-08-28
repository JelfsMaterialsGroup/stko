"""
Utilities for periodic molecules.

"""

import logging
import numpy as np

logger = logging.getLogger(__name__)


def get_approximate_cell_size(molecule, vector_1, vector_2, vector_3):
    """
    Cell size determined from projection of atoms on cell vectors.

    Arguments
    ---------
    molecule : :class:`stk.Molecule`
        Molecule to get approximate cell size of.

    vector_1 : :class:`numpy.ndarray`
        Cell lattice vector of shape (3, ) in x direction in
        Angstrom.

    vector_2 : :class:`numpy.ndarray`
        Cell lattice vector of shape (3, ) in y direction in
        Angstrom.

    vector_3 : :class:`numpy.ndarray`
        Cell lattice vector of shape (3, ) in z direction in
        Angstrom.

    Returns
    -------
    max_extent : :class:`float`
        Get maximum multiplier of cell vectors required to fit
        molecule.

    """

    def get_projection(start, target):
        """
        Get the projection of `start` onto `target`.

        """

        return target * np.dot(
            start,
            target,
        ) / np.dot(target, target)

    def get_magnitude_along_vector(atom_pos, vector):
        """
        Get magnitude of atom position along vector.

        """

        projection = get_projection(
            start=atom_pos,
            target=vector
        )

        return np.linalg.norm(np.divide(
            projection,
            vector,
            out=np.zeros_like(projection),
            where=vector != 0
        ))

    # Get projection of all atom positions on three vectors.
    extent_of_vector_1 = 0
    extent_of_vector_2 = 0
    extent_of_vector_3 = 0
    for atom_pos in molecule.get_position_matrix():
        v1_mag = get_magnitude_along_vector(
            atom_pos=atom_pos,
            vector=vector_1
        )
        v2_mag = get_magnitude_along_vector(
            atom_pos=atom_pos,
            vector=vector_2
        )
        v3_mag = get_magnitude_along_vector(
            atom_pos=atom_pos,
            vector=vector_3
        )
        extent_of_vector_1 = max([v1_mag, extent_of_vector_1])
        extent_of_vector_2 = max([v2_mag, extent_of_vector_2])
        extent_of_vector_3 = max([v3_mag, extent_of_vector_3])

    max_extent = max([extent_of_vector_1, extent_of_vector_2])*1.35

    return max_extent
