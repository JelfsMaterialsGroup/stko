import numpy as np
from scipy.spatial.distance import euclidean


def is_valid_xtb_solvent(gfn_version, solvent_model, solvent):
    """
    Check if solvent is valid for the given GFN version.

    Parameters
    ----------
    gfn_version : :class:`int`
        GFN parameterization version. Can be: ``0``, ``1`` or ``2``.

    solvent_model : class:`str`
        Solvent model being used [1]_.

    solvent : :class:`str`
        Solvent being tested [1]_.

    Returns
    -------
    :class:`bool`
        ``True`` if solvent is valid.

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    """

    if gfn_version == 0:
        return False
    elif gfn_version == 1 and solvent_model == "gbsa":
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
    elif gfn_version == 2 and solvent_model == "gbsa":
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
    elif gfn_version == 2 and solvent_model == "alpb":
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


def get_atom_distance(position_matrix, atom1_id, atom2_id):
    """
    Return the distance between two atoms.

    """

    distance = euclidean(
        u=position_matrix[atom1_id],
        v=position_matrix[atom2_id],
    )

    return float(distance)


def calculate_dihedral(pt1, pt2, pt3, pt4):
    """
    Calculate the dihedral between four points in degrees.

    Uses Praxeolitic formula --> 1 sqrt, 1 cross product
    Output in range (-180 to 180).

    From: https://stackoverflow.com/questions/20305272/
    dihedral-torsion-angle-from-four-points-in-cartesian-
    coordinates-in-python
    (new_dihedral(p))

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


def vector_angle(vector1, vector2):
    """
    Returns the angle between two vectors in radians.
    Parameters
    ----------
    vector1 : :class:`numpy.ndarray`
        The first vector.
    vector2 : :class:`numpy.ndarray`
        The second vector.
    Returns
    -------
    :class:`float`
        The angle between `vector1` and `vector2` in radians.
    """

    if np.all(np.equal(vector1, vector2)):
        return 0.0

    numerator = np.dot(vector1, vector2)
    denominator = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    # This if statement prevents returns of NaN due to floating point
    # inaccuracy.
    term = numerator / denominator
    if term >= 1.0:
        return 0.0
    if term <= -1.0:
        return np.pi
    return np.arccos(term)


def calculate_angle(pt1, pt2, pt3):
    """
    Calculate the angle between three points in degrees.

    """

    v1 = pt1 - pt2
    v2 = pt3 - pt2
    return np.degrees(vector_angle(v1, v2))


def get_torsion_info_angles(mol, torsion_info):
    """
    Get the angles for torsion_info in mol.

    The first angle returned is torsion angle in the
    :class:`stk.ConstructedMolecule`.
    The second angle returned is torsion angle in the
    :class:`stk.BuildingBlock`.
    Both angles are in degrees.

    A :class:`stko.MatchedTorsionCalculator` should yield torsions
    such that the two angles returned are the same.

    Parameters
    ----------
    mol : :class:`.ConstructedMolecule`
        The :class:`.ConstructedMolecule` for which angles are
        computed.

    torsion_info : TorsionInfo
        Specifies the torsion for which angles will be computed.

    Returns
    -------
    angle : :class:`float`, bb_angle : :class:`float`

    """

    torsion = torsion_info.get_torsion()
    angle = calculate_dihedral(
        pt1=tuple(mol.get_atomic_positions(torsion.get_atom_ids()[0]))[0],
        pt2=tuple(mol.get_atomic_positions(torsion.get_atom_ids()[1]))[0],
        pt3=tuple(mol.get_atomic_positions(torsion.get_atom_ids()[2]))[0],
        pt4=tuple(mol.get_atomic_positions(torsion.get_atom_ids()[3]))[0],
    )
    bb_torsion = torsion_info.get_building_block_torsion()
    if bb_torsion is None:
        bb_angle = None
    else:
        bb_angle = calculate_dihedral(
            pt1=tuple(
                torsion_info.get_building_block().get_atomic_positions(
                    bb_torsion.get_atom_ids()[0]
                )
            )[0],
            pt2=tuple(
                torsion_info.get_building_block().get_atomic_positions(
                    bb_torsion.get_atom_ids()[1]
                )
            )[0],
            pt3=tuple(
                torsion_info.get_building_block().get_atomic_positions(
                    bb_torsion.get_atom_ids()[2]
                )
            )[0],
            pt4=tuple(
                torsion_info.get_building_block().get_atomic_positions(
                    bb_torsion.get_atom_ids()[3]
                )
            )[0],
        )
    return angle, bb_angle
