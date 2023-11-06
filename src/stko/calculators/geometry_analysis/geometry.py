import logging

import numpy as np
import stk

logger = logging.getLogger(__name__)


class GeometryAnalyser:
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
