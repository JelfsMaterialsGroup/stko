import numpy as np
import stk
from scipy.spatial.distance import euclidean

from stko._internal.molecular.torsion.torsion_info import TorsionInfo


def is_valid_xtb_solvent(
    gfn_version: int,
    solvent_model: str,
    solvent: str,
) -> bool:
    """Check if solvent is valid for the given GFN version.

    See Also:
        * Valid solvents: https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    Parameters:
        gfn_version:
            GFN parameterization version. Can be: ``0``, ``1`` or ``2``.

        solvent_model:
            Solvent model being used.

        solvent:
            Solvent being tested.

    Returns:
        ``True`` if solvent is valid.

    """
    if gfn_version == 0:
        return False
    if gfn_version == 1 and solvent_model == "gbsa":
        valid_solvents = {
            "acetone",
            "acetonitrile",
            "benzene",
            "CH2Cl2".lower(),
            "CHCl3".lower(),
            "CS2".lower(),
            "DMSO".lower(),
            "ether",
            "H2O".lower(),
            "methanol",
            "THF".lower(),
            "toluene",
            "water",
        }
    elif gfn_version == 1 and solvent_model == "alpb":
        valid_solvents = {
            "acetone",
            "acetonitrile",
            "aniline",
            "benzaldehyde",
            "benzene",
            "CH2Cl2".lower(),
            "CHCl3".lower(),
            "CS2".lower(),
            "dioxane",
            "DMF".lower(),
            "DMSO".lower(),
            "ether",
            "ethylacetate",
            "furane",
            "hexandecane",
            "hexane",
            "H2O".lower(),
            "nitromethane",
            "octanol",
            "octanol (wet)",
            "phenol",
            "THF".lower(),
            "toluene",
            "water",
        }
    elif gfn_version == 2 and solvent_model == "gbsa":  # noqa: PLR2004
        valid_solvents = {
            "acetone",
            "acetonitrile",
            "benzene",
            "CH2Cl2".lower(),
            "CHCl3".lower(),
            "CS2".lower(),
            "DMSO".lower(),
            "ether",
            "hexane",
            "methanol",
            "H2O".lower(),
            "THF".lower(),
            "toluene",
            "water",
        }
    elif gfn_version == 2 and solvent_model == "alpb":  # noqa: PLR2004
        valid_solvents = {
            "acetone",
            "acetonitrile",
            "aniline",
            "benzaldehyde",
            "benzene",
            "CH2Cl2".lower(),
            "CHCl3".lower(),
            "CS2".lower(),
            "dioxane",
            "DMF".lower(),
            "DMSO".lower(),
            "ether",
            "ethylacetate",
            "furane",
            "hexandecane",
            "hexane",
            "H2O".lower(),
            "nitromethane",
            "octanol",
            "octanol (wet)",
            "phenol",
            "THF".lower(),
            "toluene",
            "water",
        }
    return solvent in valid_solvents


def get_atom_distance(
    position_matrix: np.ndarray,
    atom1_id: int,
    atom2_id: int,
) -> float:
    """Return the distance between two atoms."""
    distance = euclidean(
        u=position_matrix[atom1_id],
        v=position_matrix[atom2_id],
    )

    return float(distance)


def calculate_dihedral(
    pt1: np.ndarray,
    pt2: np.ndarray,
    pt3: np.ndarray,
    pt4: np.ndarray,
) -> float:
    """Calculate the dihedral between four points in degrees.

    Uses Praxeolitic formula --> 1 sqrt, 1 cross product
    Output in range (-180 to 180).

    From: https://stackoverflow.com/a/34245697

    """
    p0 = np.asarray(pt1)
    p1 = np.asarray(pt2)
    p2 = np.asarray(pt3)
    p3 = np.asarray(pt4)

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def unit_vector(vector: np.ndarray) -> np.ndarray:
    """Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)


def vector_angle(vector1: np.ndarray, vector2: np.ndarray) -> float:
    """Returns the angle between two vectors in radians.

    From: https://stackoverflow.com/a/13849249

    Parameters:
        vector1:
            The first vector.

        vector2:
            The second vector.

    Returns:
        The angle between `vector1` and `vector2` in radians.

    """
    v1_u = unit_vector(vector1)
    v2_u = unit_vector(vector2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def calculate_angle(
    pt1: np.ndarray,
    pt2: np.ndarray,
    pt3: np.ndarray,
) -> float:
    """Calculate the angle between three points in degrees."""
    v1 = pt1 - pt2
    v2 = pt3 - pt2
    return np.degrees(vector_angle(v1, v2))


def get_torsion_info_angles(
    mol: stk.ConstructedMolecule,
    torsion_info: TorsionInfo,
) -> tuple[float, float | None]:
    """Get the angles for torsion_info in mol.

    The first angle returned is torsion angle in the
    :class:`stk.ConstructedMolecule`.
    The second angle returned is torsion angle in the
    :class:`stk.BuildingBlock`.

    A :class:`stko.MatchedTorsionCalculator` should yield torsions
    such that the two angles returned are the same.

    Parameters:
        mol:
            The molecule for which angles are computed.

        torsion_info:
            Specifies the torsion for which angles will be computed.

    Returns:
        The angle and the bb_angle in degrees.

    """
    torsion = torsion_info.get_torsion()
    angle = calculate_dihedral(
        pt1=next(mol.get_atomic_positions(torsion.get_atom_ids()[0])),
        pt2=next(mol.get_atomic_positions(torsion.get_atom_ids()[1])),
        pt3=next(mol.get_atomic_positions(torsion.get_atom_ids()[2])),
        pt4=next(mol.get_atomic_positions(torsion.get_atom_ids()[3])),
    )
    bb_torsion = torsion_info.get_building_block_torsion()
    bb = torsion_info.get_building_block()
    if bb_torsion is None or bb is None:
        bb_angle = None
    else:
        bb_angle = calculate_dihedral(
            pt1=next(bb.get_atomic_positions(bb_torsion.get_atom_ids()[0])),
            pt2=next(bb.get_atomic_positions(bb_torsion.get_atom_ids()[1])),
            pt3=next(bb.get_atomic_positions(bb_torsion.get_atom_ids()[2])),
            pt4=next(bb.get_atomic_positions(bb_torsion.get_atom_ids()[3])),
        )
    return angle, bb_angle
