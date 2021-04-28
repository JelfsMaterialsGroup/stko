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
    ):

        flat_tuple = []
        for i in atom_ids_to_disconnect:
            flat_tuple.append(i[0])
            flat_tuple.append(i[1])
        flat_tuple = tuple(set(flat_tuple))

        connected_graphs = self.get_connected_graphs(
            molecule=molecule,
            atom_ids_to_disconnect=atom_ids_to_disconnect,
        )

        # Get centroids.
        self._centroids = {}
        self._connectivities = {}
        for i, cg in enumerate(connected_graphs):
            self._centroids[i] = molecule.get_centroid(
                atom_ids=[i.get_id() for i in list(cg)]
            )
            disconnections = 0
            for atom in list(cg):
                aid = atom.get_id()
                if aid in flat_tuple:
                    disconnections += 1
            self._connectivities[i] = disconnections

    def get_connected_graphs(
        self,
        molecule,
        atom_ids_to_disconnect
    ):
        graph = Network(molecule)
        graph.delete_bonds(atom_ids_to_disconnect)
        connected_graphs = graph.get_connected_graphs()
        return connected_graphs

    def get_vertex_positions(self):
        return self._centroids

    def get_connectivities(self):
        return self._connectivities

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'<{self.__class__.__name__} at {id(self)}>'
