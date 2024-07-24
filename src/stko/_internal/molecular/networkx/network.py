import logging
from collections import abc
from typing import Self

import networkx as nx
import stk

from stko._internal.molecular.atoms.positioned_atom import PositionedAtom

logger = logging.getLogger(__name__)


class Network:
    """Definition of a :mod:`networkx` graph of an :class:`stk.Molecule`.

    Parameters:
        graph:
            The NetworkX graph to initialise from.

    Examples:
        An stk molecule can be converted into a NetworkX object. This allows
        for the disconnection and manipulation of the molecular graph.

        .. code-block:: python

            import stk
            import stko

            molecule = stk.BuildingBlock('NCCNCCN').with_centroid(
                position=np.array((10, 10, 10))
            )
            graph = stko.Network.init_from_molecule(molecule)

            # Pick some bonds to break based on atom ids in those bonds.
            atom_ids_to_disconnect = ((2, 3),)
            graph = graph.with_deleted_bonds(atom_ids_to_disconnect)

            # A series of graphs still connected.
            # This gives a list of atoms in each graph.
            connected_graphs = graph.get_connected_components()

            # Get centroids of connected graphs.
            centroids = [
                molecule.get_centroid(
                    atom_ids=tuple(i.get_id() for i in list(cg))
                )
                for cg in connected_graphs
            ]

            # Get individual stk.Molecule classes of each connected graph.
            # This is what stko.molecule_analysis.DecomposeMOC does
            # after deleting metal atoms!
            for atom_ids in connected_graphs:
                # Get atoms from nodes.
                atoms = list(cg)
                atom_ids = tuple(i.get_id() for i in atoms)

                # Sort both by atom id.
                atom_ids, atoms = zip(
                    *sorted(zip(atom_ids, atoms, strict=True)), strict=True
                )

                # Map old ids to a new molecule ids.
                atom_ids_map = {atom_ids[i]: i for i in range(len(atom_ids))}

                # Make a new stk.Building Block.
                new_mol = stk.BuildingBlock.init(
                    atoms=(
                        stk.Atom(
                            id=atom_ids_map[i.get_id()],
                            atomic_number=i.get_atomic_number(),
                            charge=i.get_charge(),
                        )
                        for i in atoms
                    ),
                    bonds=(
                        i.with_ids(id_map=atom_ids_map)
                        for i in molecule.get_bonds()
                        if i.get_atom1().get_id() in atom_ids
                        and i.get_atom2().get_id() in atom_ids
                    ),
                    position_matrix=np.array(
                        tuple(
                            i
                            for i in molecule.get_atomic_positions(
                                atom_ids=atom_ids
                            )
                        )
                    ),
                )

                # Now do something with this molecule!

    """

    def __init__(self, graph: nx.Graph) -> None:
        self._graph = graph

    @classmethod
    def init_from_molecule(cls, molecule: stk.Molecule) -> Self:
        """Initialize from an stk molecule.

        Parameters:
            molecule:
                The molecule to initialise from.

        """
        g = nx.Graph()

        pos_mat = molecule.get_position_matrix()
        for atom in molecule.get_atoms():
            pos = tuple(float(i) for i in pos_mat[atom.get_id()])
            pa = PositionedAtom(atom, pos)
            g.add_node(pa)

        # Define edges.
        for bond in molecule.get_bonds():
            n1, n2 = (
                i
                for i in g.nodes
                if i.get_id()
                in (
                    bond.get_atom1().get_id(),
                    bond.get_atom2().get_id(),
                )
            )

            g.add_edge(
                n1,
                n2,
                order=bond.get_order(),
                periodicity=bond.get_periodicity(),
            )

        return cls(g)

    def get_graph(self) -> nx.Graph:
        """Return a :class:`networkx.Graph`."""
        return self._graph

    def get_nodes(self) -> abc.Iterator[PositionedAtom]:
        """Yield nodes of :class:`networkx.Graph` (:class:`PositionAtom`)."""
        yield from self._graph.nodes

    def clone(self) -> Self:
        """Return a clone."""
        clone = self.__class__.__new__(self.__class__)
        Network.__init__(self=clone, graph=self._graph)
        return clone

    def _with_deleted_bonds(
        self,
        atom_ids: abc.Iterable[tuple[int, int]],
    ) -> Self:
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

    def with_deleted_bonds(
        self,
        atom_ids: abc.Iterable[tuple[int, int]],
    ) -> Self:
        """Return a clone with edges between `atom_ids` deleted."""
        return self.clone()._with_deleted_bonds(atom_ids)  # noqa: SLF001

    def _with_deleted_elements(self, atomic_numbers: tuple[int]) -> Self:
        to_delete = []
        deleted_atom_ids = set()
        for node in self._graph.nodes:
            if node.get_atomic_number() in atomic_numbers:
                deleted_atom_ids.add(node.get_id())
                to_delete.append(node)

        self._graph.remove_nodes_from(to_delete)

        # Remove associated edges.
        to_delete = []
        for edge in self._graph.edges:
            a1id = edge[0].get_id()
            a2id = edge[1].get_id()
            if a1id in deleted_atom_ids or a2id in deleted_atom_ids:
                to_delete.append(edge)

        for id1, id2 in to_delete:
            self._graph.remove_edge(id1, id2)

        return self

    def with_deleted_elements(self, atomic_numbers: tuple[int]) -> Self:
        """Return a clone with nodes with `atomic numbers` deleted.

        .. warning::
            This code is only present in the latest versions of stko
            that require Python 3.11!

        """
        return self.clone()._with_deleted_elements(atomic_numbers)  # noqa: SLF001

    def get_connected_components(self) -> list[nx.Graph]:
        """Get connected components within full graph.

        Returns:
            List of connected components of graph.

        """
        return [
            self._graph.subgraph(c).copy()
            for c in sorted(nx.connected_components(self._graph))
        ]

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"n={self._graph.number_of_nodes()}, "
            f"e={self._graph.number_of_edges()})"
        )
