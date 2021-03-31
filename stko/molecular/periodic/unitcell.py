"""
Unit Cell
=========

Class holding periodic cell information.

"""

import numpy as np
import logging
from stk import PeriodicInfo

from .utilities import get_from_parameters

logger = logging.getLogger(__name__)


class UnitCell(PeriodicInfo):
    """
    Unit cell information for periodic systems.

    """

    @classmethod
    def _update_periodic_info(cls, vector_1, vector_2, vector_3):
        """
        Return clone of :class:`.UnitCell` with new parameters.

        """

        clone = cls.__new__(cls)
        UnitCell.__init__(
            self=clone,
            vector_1=vector_1,
            vector_2=vector_2,
            vector_3=vector_3,
        )

        return clone

    def with_cell_from_vectors(self, vector_1, vector_2, vector_3):
        """
        Update cell.

        Parameters
        ----------
        vector_1 : :class:`numpy.ndarray`
            First cell lattice vector of shape (3, ) in
            Angstrom.

        vector_2 : :class:`numpy.ndarray`
            Second cell lattice vector of shape (3, ) in
            Angstrom.

        vector_3 : :class:`numpy.ndarray`
            Third cell lattice vector of shape (3, ) in
            Angstrom.

        Returns
        -------
        :class:`.UnitCell`
            Clone with updated cell parameters.

        """

        return self.__class__._update_periodic_info(
            vector_1=vector_1,
            vector_2=vector_2,
            vector_3=vector_3,
        )

    def with_cell_from_turbomole(self, filename):
        """
        Update cell from structure in Turbomole coord file.

        Returns
        -------
        :class:`.UnitCell`
            Clone with updated cell parameters.

        """

        bohr_to_ang = 0.5291772105638411

        with open(filename, 'r') as f:
            content = f.readlines()

        periodicity = False
        lattice_vectors = None
        lattice_units = None
        cell_parameters = None
        cell_units = None
        for line_number, line in enumerate(content):
            if '$periodic' in line:
                periodicity = int(line.rstrip().split()[1])
            if '$cell' in line:
                if 'angs' in line:
                    cell_units = 'angstrom'
                elif 'bohr' in line:
                    cell_units = 'bohr'
                else:
                    raise ValueError('cell not in Angstroms.')
                cell_parameters = [
                    float(j)
                    for j in content[line_number+1].rstrip().split()
                ]
            if '$lattice' in line:
                if 'angs' in line:
                    lattice_units = 'angstrom'
                elif 'bohr' in line:
                    lattice_units = 'bohr'
                else:
                    raise ValueError('lattice not in Angstroms.')
                lattice_vectors = (
                    np.array([
                        float(j) for j
                        in content[line_number+1].rstrip().split()
                    ]),
                    np.array([
                        float(j) for j
                        in content[line_number+2].rstrip().split()
                    ]),
                    np.array([
                        float(j) for j
                        in content[line_number+3].rstrip().split()
                    ]),
                )

        # Check that cell is only defined once.
        chk2 = (
            lattice_vectors is not None and cell_parameters is not None
        )
        if periodicity and chk2:
            raise RuntimeError(
                'The cell is defined twice in the file.'
            )

        if lattice_vectors is not None:
            vector_1 = (
                lattice_vectors[0]*bohr_to_ang
                if lattice_units == 'bohr' else lattice_vectors[0]
            )
            vector_2 = (
                lattice_vectors[1]*bohr_to_ang
                if lattice_units == 'bohr' else lattice_vectors[0]
            )
            vector_3 = (
                lattice_vectors[2]*bohr_to_ang
                if lattice_units == 'bohr' else lattice_vectors[0]
            )
        elif cell_parameters is not None:
            vector_1, vector_2, vector_3 = get_from_parameters(
                a=(
                    cell_parameters[0] * bohr_to_ang
                    if cell_units == 'bohr'
                    else cell_parameters[0]
                ),
                b=(
                    cell_parameters[1] * bohr_to_ang
                    if cell_units == 'bohr'
                    else cell_parameters[1]
                ),
                c=(
                    cell_parameters[2] * bohr_to_ang
                    if cell_units == 'bohr'
                    else cell_parameters[2]
                ),
                alpha=cell_parameters[3],
                beta=cell_parameters[4],
                gamma=cell_parameters[5],
            )
        else:
            raise RuntimeError(
                'The cell is not defined in the file.'
            )
        # Update the cell.
        return self._update_periodic_info(vector_1, vector_2, vector_3)

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

        vector_1, vector_2, vector_3 = get_from_parameters(
            a=cell_info['a'],
            b=cell_info['b'],
            c=cell_info['c'],
            alpha=cell_info['alpha'],
            beta=cell_info['beta'],
            gamma=cell_info['gamma'],
        )
        # Update the cell.
        return self._update_periodic_info(vector_1, vector_2, vector_3)
