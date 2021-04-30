"""
Topology Info
=============

Class containing extracted topology information.

"""

import logging

logger = logging.getLogger(__name__)


class TopologyInfo:
    """
    Extracted information of a topology.

    """

    def __init__(self, centroids, connectivities, edge_pairs):
        """


        """

        self._centroids = centroids
        self._connectivities = connectivities
