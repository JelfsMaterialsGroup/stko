import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class TopologyInfo:
    """Extracted information of a topology.

    Parameters:
        centroids:
            Positions of vertices.

        connectivities:
            Connections between vertices.

        edge_pairs:
            Pairs of vertices with edges between them.

    """

    def __init__(
        self,
        centroids: dict,
        connectivities: dict,
        edge_pairs: list[tuple[int, int]],
    ) -> None:
        self._centroids = centroids
        self._connectivities = connectivities
        self._edge_pairs = edge_pairs

    def get_vertex_positions(self) -> dict:
        """Get the positions of each vertex.

        Returns:
            Vertex ids with their positions.

        """
        return self._centroids

    def get_connectivities(self) -> dict:
        """Get the number of connections of each vertex.

        Returns:
            Vertex ids with their number of connections.

        """
        return self._connectivities

    def get_edge_pairs(self) -> list[tuple[int, int]]:
        """Get the edge pairs.

        Returns:
            List of edge pairs.

        """
        return self._edge_pairs

    def write(self, path: Path) -> None:
        """Writes a mock .pdb with vertex centroids and edges as bonds."""
        content = []

        atom_counts: dict[str, int] = {}
        hetatm = "HETATM"
        alt_loc = ""
        res_name = "UNL"
        chain_id = ""
        res_seq = "1"
        i_code = ""
        occupancy = "1.00"
        temp_factor = "0.00"

        # This set will be used by bonds.
        atoms = set()
        for cent in self._centroids:
            atoms.add(cent)
            serial = cent + 1
            element = "P"
            charge = 0
            atom_counts[element] = atom_counts.get(element, 0) + 1
            name = f"{element}{atom_counts[element]}"
            # Make sure the coords are no more than 8 columns wide
            # each.
            x, y, z = self._centroids[cent]

            content.append(
                f"{hetatm:<6}{serial:>5} {name:<4}"
                f"{alt_loc:<1}{res_name:<3} {chain_id:<1}"
                f"{res_seq:>4}{i_code:<1}   "
                f" {x:>7.3f} {y:>7.3f} {z:>7.3f}"
                f"{occupancy:>6}{temp_factor:>6}          "
                f"{element:>2}{charge:>2}\n"
            )

        conect = "CONECT"
        for edge in self._edge_pairs:
            a1 = edge[0]
            a2 = edge[1]
            if a1 in atoms and a2 in atoms:
                content.append(
                    f"{conect:<6}{a1+1:>5}{a2+1:>5}               \n"
                )

        content.append("END\n")

        path.write_text("".join(content))
