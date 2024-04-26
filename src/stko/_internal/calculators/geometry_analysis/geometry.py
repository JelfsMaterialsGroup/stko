import itertools as it
import logging
from collections import defaultdict

import numpy as np
import stk
from rdkit.Chem import AllChem as rdkit  # noqa: N813
from scipy.spatial.distance import cdist, pdist

from stko._internal.utilities.utilities import (
    calculate_dihedral,
    get_atom_distance,
    vector_angle,
)

logger = logging.getLogger(__name__)


class GeometryAnalyser:
    """Tools for analysing the geometry of molecules.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def _get_metal_atom_ids(
        self,
        molecule: stk.Molecule,
        metal_atom_nos: tuple[int, ...],
    ) -> list[int]:
        return [
            i.get_id()
            for i in molecule.get_atoms()
            if i.get_atomic_number() in metal_atom_nos
        ]

    def get_metal_distances(
        self,
        molecule: stk.Molecule,
        metal_atom_nos: tuple[int, ...],
    ) -> dict[tuple[int, int], float]:
        """Get all metal atom pair distances.

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
        metal_atom_nos: tuple[int, ...],
    ) -> dict[tuple[int, int], float]:
        """Get all metal-centroid-metal angles.

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
        """Get the minimum distance between the molecule and centroid.

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
        """Get the average distance between the molecule and centroid.

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

    def _get_paths(
        self,
        molecule: stk.Molecule,
        path_length: int,
    ) -> tuple[tuple[int, ...], ...]:
        return rdkit.FindAllPathsOfLengthN(
            mol=molecule.to_rdkit_mol(),
            length=path_length,
            useBonds=False,
            useHs=True,
        )

    def get_min_atom_atom_distance(self, molecule: stk.Molecule) -> float:
        """Get the minimum distance between atoms in molecule.

        This does not consider bonding.

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            The minimum distance.

        """
        pair_dists = pdist(molecule.get_position_matrix())
        return np.min(pair_dists.flatten())

    def get_radius_gyration(self, molecule: stk.Molecule) -> float:
        """Get the radius of gyration of the molecule.

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            R_g in Angstrom.

        """
        centroid = molecule.get_centroid()
        pos_mat = molecule.get_position_matrix()
        vectors = pos_mat - centroid
        distances2 = np.square(np.linalg.norm(vectors, axis=1))
        rg2 = (1 / molecule.get_num_atoms()) * np.sum(distances2)
        return np.sqrt(rg2)

    def get_max_diameter(self, molecule: stk.Molecule) -> float:
        """Get the maximum diameter of the molecule (defined in stk).

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            The maximum diameter in Angstrom.

        """
        return molecule.get_maximum_diameter()

    def calculate_bonds(
        self,
        molecule: stk.Molecule,
    ) -> dict[tuple[str, str], list[float]]:
        """Calculate bond lengths for all `stk.Molecule.get_bonds()`.

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            Dictionary of bonds organised by element pair.

        """
        position_matrix = molecule.get_position_matrix()
        lengths: dict[tuple[str, str], list[float]] = defaultdict(list)
        for bond in molecule.get_bonds():
            a1id = bond.get_atom1().get_id()
            a2id = bond.get_atom2().get_id()
            a, b = sorted(
                (
                    bond.get_atom1().__class__.__name__,
                    bond.get_atom2().__class__.__name__,
                )
            )
            lengths[(a, b)].append(
                get_atom_distance(position_matrix, a1id, a2id)
            )

        return lengths

    def calculate_angles(
        self,
        molecule: stk.Molecule,
    ) -> dict[tuple[str, str, str], list[float]]:
        """Calculate angles for all angles defined by molecule bonding.

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            Dictionary of angles organised by element triplet.

        """
        position_matrix = molecule.get_position_matrix()
        angles: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        for a_ids in self._get_paths(molecule, 3):
            atoms = list(molecule.get_atoms(atom_ids=a_ids))
            atom1 = atoms[0]
            atom2 = atoms[1]
            atom3 = atoms[2]
            angle_type_option1 = (
                atom1.__class__.__name__,
                atom2.__class__.__name__,
                atom3.__class__.__name__,
            )
            angle_type_option2 = (
                atom3.__class__.__name__,
                atom2.__class__.__name__,
                atom1.__class__.__name__,
            )

            vector1 = (
                position_matrix[atom2.get_id()]
                - position_matrix[atom1.get_id()]
            )
            vector2 = (
                position_matrix[atom2.get_id()]
                - position_matrix[atom3.get_id()]
            )

            if angle_type_option1 in angles:
                angles[angle_type_option1].append(
                    np.degrees(vector_angle(vector1, vector2))
                )
            elif angle_type_option2 in angles:
                angles[angle_type_option2].append(
                    np.degrees(vector_angle(vector1, vector2))
                )
            else:
                angles[angle_type_option1].append(
                    np.degrees(vector_angle(vector1, vector2))
                )

        return angles

    def calculate_torsions(
        self,
        molecule: stk.Molecule,
    ) -> dict[tuple[str, ...], list[float]]:
        """Calculate all torsions defined by molecule bonding.

        Parameters:
            molecule:
                The molecule to analyse.

        Returns:
            Dictionary of torsions organised by elements.

        """
        position_matrix = molecule.get_position_matrix()

        torsions: dict[tuple[str, ...], list[float]] = defaultdict(list)
        for a_ids in self._get_paths(molecule, 4):
            atoms = list(molecule.get_atoms(atom_ids=a_ids))
            atom1 = atoms[0]
            atom2 = atoms[1]
            atom3 = atoms[2]
            atom4 = atoms[3]
            torsion_type_option1 = (
                atom1.__class__.__name__,
                atom2.__class__.__name__,
                atom3.__class__.__name__,
                atom4.__class__.__name__,
            )
            torsion_type_option2 = (
                atom4.__class__.__name__,
                atom3.__class__.__name__,
                atom2.__class__.__name__,
                atom1.__class__.__name__,
            )

            if torsion_type_option1 in torsions:
                torsions[torsion_type_option1].append(
                    calculate_dihedral(
                        pt1=position_matrix[atom1.get_id()],
                        pt2=position_matrix[atom2.get_id()],
                        pt3=position_matrix[atom3.get_id()],
                        pt4=position_matrix[atom4.get_id()],
                    )
                )
            elif torsion_type_option2 in torsions:
                torsions[torsion_type_option2].append(
                    calculate_dihedral(
                        pt1=position_matrix[atom4.get_id()],
                        pt2=position_matrix[atom3.get_id()],
                        pt3=position_matrix[atom2.get_id()],
                        pt4=position_matrix[atom1.get_id()],
                    )
                )
            else:
                torsions[torsion_type_option1].append(
                    calculate_dihedral(
                        pt1=position_matrix[atom1.get_id()],
                        pt2=position_matrix[atom2.get_id()],
                        pt3=position_matrix[atom3.get_id()],
                        pt4=position_matrix[atom4.get_id()],
                    )
                )

        return torsions
