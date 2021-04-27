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

    def __init__(self, vertex_options):
        self._vertex_options = vertex_options

    def extract_topology(self, molecule, bond_ids_to_break):
        raise NotImplementedError()

    def get_disconnected_graphs(self):
        graph = Network(molecule)
        graph.delete_bonds(bonds_to_break)
        return disconnected_graphs

    def get_vertices(self):
        raise NotImplementedError()

    def get_edges(self):
        raise NotImplementedError()

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'<{self.__class__.__name__} at {id(self)}>'
