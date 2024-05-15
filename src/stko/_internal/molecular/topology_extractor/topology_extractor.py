import logging
from collections import abc

import stk

from stko._internal.molecular.networkx.network import Network
from stko._internal.molecular.topology_extractor.topology_info import (
    TopologyInfo,
)

logger = logging.getLogger(__name__)


class TopologyExtractor:
    """Extractor of topology definitions from a molecule.

    Examples:
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
        molecule: stk.Molecule,
        broken_bonds_by_id: abc.Iterable[tuple[int, int]],
        disconnectors: set[int],
    ) -> TopologyInfo:
        """Extract a toplogy defining a molecule with disconnections.

        Parameters:
            molecule:
                Molecule to get underlying topology of.

            broken_bonds_by_id:
                Tuples of bonds to break by atom id.

            disconnectors:
                Atom ids of disconnection points.

        Returns:
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
                i
                for i, cg in enumerate(connected_graphs)
                if pair[0] in {atom.get_id() for atom in cg}
                or pair[1] in {atom.get_id() for atom in cg}
            )
            edge_pairs.append((a1_g, a2_g))

        return TopologyInfo(centroids, connectivities, edge_pairs)

    def get_connected_graphs(
        self,
        molecule: stk.Molecule,
        atom_ids_to_disconnect: abc.Iterable[tuple[int, int]],
    ) -> list:
        graph = Network.init_from_molecule(molecule)
        graph = graph.with_deleted_bonds(atom_ids_to_disconnect)
        return graph.get_connected_components()
