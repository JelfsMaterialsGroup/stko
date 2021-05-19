"""
Topology Info
=============

#. :class:`.TopologyInfo`

Class containing extracted topology information.

"""

import logging

logger = logging.getLogger(__name__)


class TopologyInfo:
    """
    Extracted information of a topology.

    """

    def __init__(self, centroids, connectivities, edge_pairs):
        """


        """

        self._centroids = centroids
        self._connectivities = connectivities
        self._edge_pairs = edge_pairs

    def get_vertex_positions(self):
        """
        Get the positions of each vertex.

        Returns
        -------
        :class:`dict`
            Vertex ids with their positions.

        """

        return self._centroids

    def get_connectivities(self):
        """
        Get the number of connections of each vertex.

        Returns
        -------
        :class:`dict`
            Vertex ids with their number of connections.

        """

        return self._connectivities

    def get_edge_pairs(self):
        """
        Get the edge pairs.

        Returns
        -------
        :class:`list` of :class:`tuple`
            List of edge pairs.

        """

        return self._edge_pairs

    def write(self, path):
        """
        Writes a mock .pdb with vertex centroids and edges as bonds.

        """

        content = []

        atom_counts = {}
        hetatm = 'HETATM'
        alt_loc = ''
        res_name = 'UNL'
        chain_id = ''
        res_seq = '1'
        i_code = ''
        occupancy = '1.00'
        temp_factor = '0.00'

        # This set will be used by bonds.
        atoms = set()
        for cent in self._centroids:
            atoms.add(cent)
            serial = cent+1
            element = 'P'
            charge = 0
            atom_counts[element] = atom_counts.get(element, 0) + 1
            name = f'{element}{atom_counts[element]}'
            # Make sure the coords are no more than 8 columns wide
            # each.
            x, y, z = self._centroids[cent]

            content.append(
                f'{hetatm:<6}{serial:>5} {name:<4}'
                f'{alt_loc:<1}{res_name:<3} {chain_id:<1}'
                f'{res_seq:>4}{i_code:<1}   '
                f' {x:>7.3f} {y:>7.3f} {z:>7.3f}'
                f'{occupancy:>6}{temp_factor:>6}          '
                f'{element:>2}{charge:>2}\n'
            )

        conect = 'CONECT'
        for edge in self._edge_pairs:
            a1 = edge[0]
            a2 = edge[1]
            if a1 in atoms and a2 in atoms:
                content.append(
                    f'{conect:<6}{a1+1:>5}{a2+1:>5}               \n'
                )

        content.append('END\n')

        with open(path, 'w') as f:
            f.write(''.join(content))

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'<{self.__class__.__name__} at {id(self)}>'
