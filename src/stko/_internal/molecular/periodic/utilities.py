import logging

import numpy as np
import stk

logger = logging.getLogger(__name__)


def get_approximate_cell_size(
    molecule: stk.Molecule,
    vector_1: np.ndarray,
    vector_2: np.ndarray,
    vector_3: np.ndarray,
) -> float:
    """Cell size determined from projection of atoms on cell vectors.

    Parameters:
        molecule:
            Molecule to get approximate cell size of.

        vector_1:
            First cell lattice vector of shape (3, ) in Angstrom.

        vector_2:
            Second cell lattice vector of shape (3, ) in Angstrom.

        vector_3:
            Third cell lattice vector of shape (3, ) in Angstrom.

    Returns:
        The maximum multiplier of cell vectors required to fit
        molecule.

    """

    def get_projection(start: np.ndarray, target: np.ndarray) -> np.ndarray:
        """Get the projection of `start` onto `target`."""
        return (
            target
            * np.dot(
                start,
                target,
            )
            / np.dot(target, target)
        )

    def get_magnitude_along_vector(
        atom_pos: np.ndarray,
        vector: np.ndarray,
    ) -> float:
        """Get magnitude of atom position along vector."""
        projection = get_projection(start=atom_pos, target=vector)

        return float(
            np.linalg.norm(
                np.divide(
                    projection,
                    vector,
                    out=np.zeros_like(projection),
                    where=vector != 0,
                )
            )
        )

    # Get projection of all atom positions on three vectors.
    extent_of_vector_1 = 0.0
    extent_of_vector_2 = 0.0
    extent_of_vector_3 = 0.0
    for atom_pos in molecule.get_position_matrix():
        v1_mag = get_magnitude_along_vector(atom_pos=atom_pos, vector=vector_1)
        v2_mag = get_magnitude_along_vector(atom_pos=atom_pos, vector=vector_2)
        v3_mag = get_magnitude_along_vector(atom_pos=atom_pos, vector=vector_3)
        extent_of_vector_1 = max([v1_mag, extent_of_vector_1])
        extent_of_vector_2 = max([v2_mag, extent_of_vector_2])
        extent_of_vector_3 = max([v3_mag, extent_of_vector_3])

    return max([extent_of_vector_1, extent_of_vector_2]) * 1.35


def get_from_parameters(  # noqa: PLR0913
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> tuple[np.ndarray, ...]:
    """Create a Lattice using unit cell lengths and angles (in degrees).

    This code is modified from the pymatgen source code [1]_.

    Parameters:
        a:
            *a* lattice parameter.

        b:
            *b* lattice parameter.

        c:
            *c* lattice parameter.

        alpha:
            *alpha* angle in degrees.

        beta:
            *beta* angle in degrees.

        gamma:
            *gamma* angle in degrees.


    Returns:
        Tuple of cell lattice vectors of shape (3, ) in Angstrom.

    References:
        .. [1] https://pymatgen.org/pymatgen.util.num.html

    """
    angles_r = np.radians([alpha, beta, gamma])
    cos_alpha, cos_beta, cos_gamma = np.cos(angles_r)
    sin_alpha, sin_beta, _ = np.sin(angles_r)

    val = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
    # Sometimes rounding errors result in values slightly > 1.
    val = cap_absolute_value(val)
    gamma_star = np.arccos(val)

    vector_a = np.array([a * sin_beta, 0.0, a * cos_beta])
    vector_b = np.array(
        [
            -b * sin_alpha * np.cos(gamma_star),
            b * sin_alpha * np.sin(gamma_star),
            b * cos_alpha,
        ]
    )
    vector_c = np.array([0.0, 0.0, float(c)])

    return (vector_a, vector_b, vector_c)


def cap_absolute_value(value: float, max_absolute_value: float = 1) -> float:
    """Returns `value` with absolute value capped at `max_absolute_value`.

    Particularly useful in passing values to trignometric functions
    where numerical errors may result in an argument > 1 being passed
    in.

    This code is modified from the pymatgen source code [2]_.

    Parameters:
        value:
            Value to cap.

        max_absolute_value:
            Absolute value to cap `value` at. Defaults to 1.

    Returns:
        `value` capped at `max_absolute_value` with sign preserved.

    References:
        .. [2] https://pymatgen.org/pymatgen.util.num.html

    """
    return max(min(value, max_absolute_value), -max_absolute_value)
