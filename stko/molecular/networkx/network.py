"""
Network
=======

Class for defining a networkx graph from a molecule.

"""

import logging

logger = logging.getLogger(__name__)


class Network:
    """
    Definition of a network of an stk.Molecule.

    """

    def __init__(self, molecule):
        self._molecule = molecule
        self._graph = self.get_graph()

    def get_graph(self):
        raise NotImplementedError()

    def delete_bonds(self, bonds_to_delete):
        raise NotImplementedError()

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'<{self.__class__.__name__} at {id(self)}>'
