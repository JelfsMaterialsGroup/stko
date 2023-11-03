import logging

import stk
import numpy as np


logger = logging.getLogger(__name__)


class ConstructedAnalyser:
    """
    Analyses geometry of building blocks in constructed molecules.

    """

    def get_building_block_centroids(
        self,
        molecule: stk.ConstructedMolecule,
    ) -> tuple[np.ndarray, ...]:
        """
        Get the centroids of all building blocks.

        Parameters:

            molecule:
                Molecule to analyse.

        """

        raise NotImplementedError()
