import logging
from collections import defaultdict

import numpy as np
import stk

logger = logging.getLogger(__name__)


class ConstructedAnalyser:
    """Analyses geometry of building blocks in constructed molecules.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def get_building_block_atom_ids(
        self,
        molecule: stk.ConstructedMolecule,
    ) -> dict[int | None, list[int]]:
        """Get the centroids of all building blocks.

        Parameters:
            molecule:
                Molecule to analyse.

        """
        atom_ids = defaultdict(list)
        for atom in molecule.get_atom_infos():
            if atom.get_building_block_id() is None:
                continue
            bb_id = atom.get_building_block_id()
            atom_ids[bb_id].append(atom.get_atom().get_id())
        return atom_ids

    def get_building_block_centroids(
        self,
        molecule: stk.ConstructedMolecule,
    ) -> dict[int | None, np.ndarray]:
        """Get the centroids of all building blocks.

        Parameters:
            molecule:
                Molecule to analyse.

        """
        atom_ids = self.get_building_block_atom_ids(molecule)
        return {
            i: molecule.get_centroid(atom_ids=atom_ids[i]) for i in atom_ids
        }
