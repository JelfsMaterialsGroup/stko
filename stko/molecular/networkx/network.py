"""
Network
=======

#. :class:`.Network`

Class for defining a :mod:`networkx` graph from a molecule.

"""

import logging
import networkx as nx

logger = logging.getLogger(__name__)


class PositionedAtom:
    """
    A container for stk.Atom and a coordinate.

    """

    def __init__(self, atom, position):
        """
        Initialize a :class:`PositionedAtom`.

        Parameters
        ----------
        atom : :class:`stk.Atom`
            The atom.

        position : : :class:`tuple` of :class:`float`
            The position (`x`, `y`, `z`) of the atom in cartesian
            coordinates.

        """

        self._atom = atom
        self._position = position

    def get_atom(self):
        return self._atom

    def get_atomic_number(self):
        return self._atom.get_atomic_number()

    def get_charge(self):
        return self._atom.get_charge()

    def get_id(self):
        return self._atom.get_id()

    def get_position(self):
        return self._position

    def _with_id(self, id):
        """
        Modify the atom id.

        """

        self._atom = self._atom.with_id(id)
        return self

    def with_id(self, id):
        """
        Get a clone but with a different id.

        Returns
        -------
        :class:`.Atom`
            A clone with a new id. Has the same type as the original
            atom.

        """

        return self.clone()._with_id(id)

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.Atom`
            The clone. It has the same type as the original atom.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._atom = self._atom
        clone._position = self._position
        return clone

    def __repr__(self):
        return (
            f'{self._atom.__class__.__name__}({self._atom.get_id()})'
        )

    def __str__(self):
        return repr(self)


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

        pos_mat = molecule.get_position_matrix()
        for atom in molecule.get_atoms():
            pos = tuple(float(i) for i in pos_mat[atom.get_id()])
            pa = PositionedAtom(atom, pos)
            g.add_node(pa)

        # Define edges.
        for bond in molecule.get_bonds():
            n1, n2 = [
                i for i in g.nodes
                if i.get_id() in (
                    bond.get_atom1().get_id(),
                    bond.get_atom2().get_id(),
                )
            ]

            g.add_edge(
                n1, n2,
                order=bond.get_order(),
                periodicity=bond.get_periodicity(),
            )

        return cls(g)

    def get_graph(self):
        """
        Return a networkx.graph.

        """

        return self._graph

    def get_nodes(self):
        """
        Yield nodes of networkx.graph.

        """

        for i in self._graph.nodes:
            yield i

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
            self._graph.subgraph(c).copy()
            for c in sorted(
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
