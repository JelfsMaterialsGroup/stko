"""
Network
=======

#. :class:`.Network`

Class for defining a :mod:`networkx` graph from a molecule.

"""

import logging
import networkx as nx

logger = logging.getLogger(__name__)


class Network:
    """
    Definition of a network of an stk.Molecule.

    """

    def __init__(self, graph):
        """
        Initialize a Network from a networkx.graph.

        """

        self._graph = graph

    @classmethod
    def init_from_molecule(cls, molecule):
        """
        Initialize a Network from a stk.Molecule.

        """

        g = nx.Graph()
        # Define by adding edges.
        for bond in molecule.get_bonds():
            g.add_edge(bond.get_atom1(), bond.get_atom2())

        return cls(g)

    def get_graph(self):
        """
        Return a networkx.graph.

        """

        return self._graph

    def clone(self):
        """
        Return a clone.

        """

        clone = self.__class__.__new__(self.__class__)
        Network.__init__(self=clone, graph=self._graph)
        return clone

    def _with_deleted_bonds(self, atom_ids):
        sorted_set = {tuple(sorted(i)) for i in atom_ids}
        to_delete = []
        for edge in self._graph.edges:
            a1id = edge[0].get_id()
            a2id = edge[1].get_id()
            pair = tuple(sorted((a1id, a2id)))
            if pair in sorted_set:
                to_delete.append(edge)

        for id1, id2 in to_delete:
            self._graph.remove_edge(id1, id2)

        return self

    def with_deleted_bonds(self, atom_ids):
        """
        Return a clone with edges between `atom_ids` deleted.

        """

        return self.clone()._with_deleted_bonds(atom_ids)

    def get_connected_components(self):
        """
        Get connected components within full graph.

        Returns
        -------
        :class:`list` of :class:`networkx.graph`
            List of connected components of graph.

        """

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
