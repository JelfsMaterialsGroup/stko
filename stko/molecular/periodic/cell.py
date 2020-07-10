"""
Cell
====

Class holding periodic cell information.

"""

import logging

logger = logging.getLogger(__name__)


class Cell:
    """
    Cell information for periodic systems.

    """

    def __init__(self, a, b, c, alpha, beta, gamma):
        """
        Initialize a :class:`MetalOptimizer` instance.

        Parameters
        ----------
        a : :class:`float`
            Length of `a` vector in Angstrom.

        b : :class:`float`
            Length of `b` vector in Angstrom.

        c : :class:`float`
            Length of `c` vector in Angstrom.

        alpha : :class:`float`
            Angle between `a` and `c` in Degrees.

        beta : :class:`float`
            Angle between `b` and `c` in Degrees

        gamma : :class:`float`
            Angle between `a` and `b` in Degrees

        """

        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

    def update_from_cif(self, filename):
        """
        Update cell from structure in CIF.

        Returns
        -------
        :class:`NoneType`
            Updates cell parameters.

        """

        cell_info = {}

        targets = {
            '_cell_length_a': 'a',
            '_cell_length_b': 'b',
            '_cell_length_c': 'c',
            '_cell_angle_alpha': 'alpha',
            '_cell_angle_beta': 'beta',
            '_cell_angle_gamma': 'gamma',
        }

        with open(filename, 'r') as f:
            lines = f.readlines()

        for targ in targets:
            for line in lines:
                # Avoid running through the rest.
                if targets[targ] in cell_info.keys():
                    break
                splits = line.rstrip().split(' ')
                if splits[0] == targ:
                    cell_info[targets[targ]] = float(splits[-1])

        self.a = cell_info['a']
        self.b = cell_info['b']
        self.c = cell_info['c']
        self.alpha = cell_info['alpha']
        self.beta = cell_info['beta']
        self.gamma = cell_info['gamma']

    def to_matrix(self):
        """
        Return matrix from of cell.

        """
        raise NotImplementedError()

    def update_from_matrix(self):
        """
        Update cell from a matrix.

        """
        raise NotImplementedError()

    def __str__(self):
        return (
            f'{self.__class__.__name__}(a={self.a}, b={self.b}, '
            f'c={self.c}, alpha={self.alpha}, beta={self.beta}, '
            f'gamma={self.gamma})'
        )

    def __repr__(self):
        return str(self)
