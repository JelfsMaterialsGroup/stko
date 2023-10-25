"""
Topology Extractor
==================

#. :class:`.TopologyExtractor`

Class for defining a topology from a molecule and disconnections.

"""

import logging

from ..networkx import Network
from .topology_info import TopologyInfo

logger = logging.getLogger(__name__)


class TopologyExtractor:
    """
    Extractor of topology definitions from a molecule.

    Examples
    --------

    Using a SMARTS string and the
    :class:`stk.SmartsFunctionalGroupFactory`, you can split a molecule
    at any point to define a topology.

    .. code-block:: python

        import stk
        import stko

        smarts = '[#6]~[#7]'

        mol = stk.BuildingBlock(
            smiles='C1=CC=C(C=C1)CN',
            functional_groups=(stk.SmartsFunctionalGroupFactory(
                smarts=smarts,
                bonders=(),
                deleters=(),
            ), )
        )

        broken_bonds_by_id = []
        disconnectors = []
        for fg in mol.get_functional_groups():
            print(fg)
            atom_ids = list(fg.get_atom_ids())
            bond_c = atom_ids[0]
            bond_n = atom_ids[1]
            broken_bonds_by_id.append(sorted((bond_c, bond_n)))
            disconnectors.extend((bond_c, bond_n))

        print(broken_bonds_by_id)
        print(disconnectors)
        print('--')

        new_topology_graph = stko.TopologyExtractor()
        tg_info = new_topology_graph.extract_topology(
            molecule=mol,
            broken_bonds_by_id=broken_bonds_by_id,
            disconnectors=set(disconnectors),
        )
        print(tg_info.get_vertex_positions())
        print(tg_info.get_connectivities())
        print(tg_info.get_edge_pairs())

    """

    def extract_topology(
        self,
        molecule,
        broken_bonds_by_id,
        disconnectors,
    ):
        """
        Extract a toplogy defining a molecule with disconnections.

        Parameters
        ----------
        molecule : :class:`stk.Molecule`
            Molecule to get underlying topology of.

        broken_bonds_by_id : :class:`iterable` of :class:`tuple`
            Tuples of bonds to break by atom id.

        disconnectors : :class:`set`
            Atom ids of disconnection points.

        Returns
        -------
        :class:`.TopologyInfo`
            Information of the underlying topology.

        """

        connected_graphs = self.get_connected_graphs(
            molecule=molecule,
            atom_ids_to_disconnect=broken_bonds_by_id,
        )

        centroids = {}
        connectivities = {}
        edge_pairs = []
        for i, cg in enumerate(connected_graphs):
            centroids[i] = molecule.get_centroid(
                atom_ids=[i.get_id() for i in cg]
            )
            disconnections = 0
            for atom in cg:
                if atom.get_id() in disconnectors:
                    disconnections += 1
            connectivities[i] = disconnections

        for pair in broken_bonds_by_id:
            a1_g, a2_g = (
                i for i, cg in enumerate(connected_graphs)
                if pair[0] in set(atom.get_id() for atom in cg)
                or pair[1] in set(atom.get_id() for atom in cg)
            )
            edge_pairs.append((a1_g, a2_g))

        return TopologyInfo(centroids, connectivities, edge_pairs)

    def get_connected_graphs(
        self,
        molecule,
        atom_ids_to_disconnect
    ):
        graph = Network.init_from_molecule(molecule)
        graph = graph.with_deleted_bonds(atom_ids_to_disconnect)
        connected_graphs = graph.get_connected_components()
        return connected_graphs

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'<{self.__class__.__name__} at {id(self)}>'
