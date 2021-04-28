"""
Network
=======

Class for defining a networkx graph from a molecule.

"""

import logging
import networkx as nx

logger = logging.getLogger(__name__)


class Network:
    """
    Definition of a network of an stk.Molecule.

    """

    def __init__(self, molecule):
        self._molecule = molecule
        self._graph = self.get_graph()

    def get_graph(self):
        g = nx.Graph()
        # Define by adding edges.
        for bond in self._molecule.get_bonds():
            g.add_edge(bond.get_atom1(), bond.get_atom2())

        return g

    def delete_bonds(self, atom_ids_to_disconnect):
        atom_ids_to_disconnect = [
            tuple(sorted(i)) for i in atom_ids_to_disconnect
        ]
        to_delete = []
        for edge in self._graph.edges:
            a1id = edge[0].get_id()
            a2id = edge[1].get_id()
            pair = tuple(sorted((a1id, a2id)))
            if pair in atom_ids_to_disconnect:
                to_delete.append(edge)

        for tdel in to_delete:
            self._graph.remove_edge(u=tdel[0], v=tdel[1])

    def get_connected_graphs(self):
        return [
            sorted(subgraph, key=lambda a: a.get_id())
            for subgraph in sorted(
                nx.connected_components(self._graph)
            )
        ]

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'n={self._graph.number_of_nodes()}, '
            f'e={self._graph.number_of_edges()})'
        )
