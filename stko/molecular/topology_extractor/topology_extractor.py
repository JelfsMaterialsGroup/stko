"""
Topology Extractor
==================

Class for defining a topology from a molecule and disconnections.

"""

import logging

from ..networkx import Network

logger = logging.getLogger(__name__)


class TopologyExtractor:
    """
    Extractor of topology definitions from a molecule.

    """

    def extract_topology(
        self,
        molecule,
        atom_ids_to_disconnect,
        vertex_options,
    ):
        connected_graphs = self.get_connected_graphs(
            molecule=molecule,
            atom_ids_to_disconnect=atom_ids_to_disconnect,
        )
        print(connected_graphs)

        # Get centroids.
        centroids = {
            i: molecule.get_centroid(
                atom_ids=[i.get_id() for i in list(cg)]
            )
            for i, cg in enumerate(connected_graphs)
        }
        print(centroids)

        # Define vertices and edges.
        self._vertices = {
            i: vertex_options[i](
                id=i,
                position=centroids[i],
                flip
            )
            for i in centroids
        }
        print(self._vertices)
        import sys
        sys.exit()
        self._edges = {
        }
        print(self._edges)
        import sys
        sys.exit()

    def get_connected_graphs(
        self,
        molecule,
        atom_ids_to_disconnect
    ):
        graph = Network(molecule)
        graph.delete_bonds(atom_ids_to_disconnect)
        connected_graphs = graph.get_connected_graphs()
        return connected_graphs

    def get_vertices(self):
        return self._vertices

    def get_edges(self):
        return self._edges

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'<{self.__class__.__name__} at {id(self)}>'
