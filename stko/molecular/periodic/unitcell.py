"""
Unit Cell
====

Class holding periodic cell information.

"""

import logging
from stk import PeriodicInfo
import numpy as np

logger = logging.getLogger(__name__)


class UnitCell(PeriodicInfo):
    """
    Unit cell information for periodic systems.

    """

    def _update_periodic_info(self, x_vector, y_vector, z_vector):
        """
        Return clone of :class:`.UnitCell` with new parameters.

        """

        clone = self.__class__.__new__(self.__class__)
        UnitCell.__init__(
            self=clone,
            x_vector=x_vector,
            y_vector=y_vector,
            z_vector=z_vector,
        )

        return clone

    def with_cell_from_turbomole(self, filename):
        """
        Update cell from structure in Turbomole coord file.

        Returns
        -------
        :class:`.UnitCell`
            Clone with updated cell parameters.

        """

        raise NotImplementedError()
        print(unit_cell)

    def with_cell_from_cif(self, filename):
        """
        Update cell from structure in CIF.

        Returns
        -------
        :class:`NoneType`
            Clone with updated cell parameters.

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

        raise NotImplementedError()

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
