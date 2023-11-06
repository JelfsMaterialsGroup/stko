import logging

import stk
import numpy as np
from collections import defaultdict


logger = logging.getLogger(__name__)


class ConstructedAnalyser:
    """
    Analyses geometry of building blocks in constructed molecules.

    """

    def get_building_block_atom_ids(
        self,
        molecule: stk.ConstructedMolecule,
    ) -> dict[int, list[int]]:
        """
        Get the centroids of all building blocks.

        Parameters:

            molecule:
                Molecule to analyse.

        """

        atom_ids = defaultdict(list)
        for atom in molecule.get_atom_infos():
            atom_ids[atom.get_building_block_id()].append(
                atom.get_atom().get_id()
            )
        return atom_ids

    def get_building_block_centroids(
        self,
        molecule: stk.ConstructedMolecule,
    ) -> dict[int, np.ndarray]:
        """
        Get the centroids of all building blocks.

        Parameters:

            molecule:
                Molecule to analyse.

        """

        atom_ids = self.get_building_block_atom_ids(molecule)
        centroids = {
            i: molecule.get_centroid(atom_ids=atom_ids[i]) for i in atom_ids
        }
        return centroids
