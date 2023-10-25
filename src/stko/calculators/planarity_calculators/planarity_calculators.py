"""
Planarity Calculators
=====================

#. :class:`.PlanarityCalculator`

Methods to calculate planarity measures of a molecule.

"""

import logging
import numpy as np

from ..calculators import Calculator
from ..results import PlanarityResults

logger = logging.getLogger(__name__)


class PlanarityCalculator(Calculator):
    """
    Calculates measures of planarity of a molecule.

    Measures based on plane deviation from Angew. paper [1]_ and a
    ChemRxiv paper [2]_.

    Plane deviation: sum of the shortest distance to the plane of best
    fit of all deviation atoms (sum abs(d_i)).

    Plane deviation span: d_max - d_min (SDP in [2]_)

    Planarity parameter: defined as
    sqrt((1/num_atoms) * (sum d_i ** 2)) (MPP in [2]_)

    Examples
    --------

    .. code-block:: python

        import stk
        import stko

        # Create a molecule whose torsions we want to know.
        mol1 = stk.BuildingBlock('c1ccccc1')

        # Create the calculator.
        pc = stko.PlanarityCalculator()

        # Extract the measures.
        pc_results = pc.get_results(mol1)
        plane_deviation = pc_results.get_plane_deviation()
        plane_deviation_span = pc_results.get_plane_deviation_span()
        planarity_parameter = pc_results.get_planarity_parameter()

    References
    ----------
    .. [1] https://onlinelibrary.wiley.com/doi/10.1002/anie.202106721

    .. [2] https://chemrxiv.org/engage/chemrxiv/article-details/
    60c73cbf9abda2e0c5f8b5c6

    """

    def _get_plane_of_best_fit(self, mol, plane_atom_ids):
        centroid = mol.get_centroid(atom_ids=plane_atom_ids)
        normal = mol.get_plane_normal(atom_ids=plane_atom_ids)
        # Plane of equation ax + by + cz = d.
        atom_plane = np.append(normal, np.sum(normal*centroid))
        return atom_plane

    def _shortest_distance_to_plane(self, plane, point):
        """
        Calculate the perpendicular distance from a point and a plane.

        """

        top = (
            plane[0]*point[0] + plane[1]*point[1] +
            plane[2]*point[2] - plane[3]
        )
        bottom = np.sqrt(plane[0]**2 + plane[1]**2 + plane[2]**2)
        distance = top / bottom
        return distance

    def _calculate_deviations(
        self,
        mol,
        atom_plane,
        deviation_atom_ids,
    ):
        return [
            self._shortest_distance_to_plane(
                plane=atom_plane,
                point=tuple(
                    mol.get_atomic_positions(atom_ids=i.get_id()),
                )[0],
            )
            for i in mol.get_atoms()
            if i.get_id() in deviation_atom_ids
        ]

    def _plane_deviation(self, deviations):
        deviations = [abs(i) for i in deviations]
        return sum(deviations)

    def _plane_deviation_span(self, deviations):
        return max(deviations) - min(deviations)

    def _planarity_parameter(self, deviations):
        deviations = [abs(i) for i in deviations]
        num_atoms = len(deviations)
        inv_num_atoms = 1/num_atoms
        sum_squared = sum(i**2 for i in deviations)
        return np.sqrt(inv_num_atoms * sum_squared)

    def calculate(
        self,
        mol,
        plane_atom_ids=None,
        deviation_atom_ids=None,
    ):
        """
        Perform calculation on `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose planarity is to be calculated.

        plane_atom_ids : iterable of :class:`int`, optional
            The atom ids to use to define the plane of best fit.

        deviation_atom_ids : iterable of :class:`int`, optional
            The atom ids to use to calculate planarity.

        Yields
        ------
        :class:`function`
            The function to perform the calculation.

        """

        if plane_atom_ids is None:
            plane_atom_ids = list(range(len(list(mol.get_atoms()))))
        else:
            plane_atom_ids = list(plane_atom_ids)

        if deviation_atom_ids is None:
            deviation_atom_ids = list(
                range(len(list(mol.get_atoms())))
            )
        else:
            deviation_atom_ids = list(deviation_atom_ids)

        atom_plane = self._get_plane_of_best_fit(mol, plane_atom_ids)
        deviations = self._calculate_deviations(
            mol=mol,
            atom_plane=atom_plane,
            deviation_atom_ids=deviation_atom_ids,
        )
        yield {
            'plane_deviation': (
                self._plane_deviation(deviations)
            ),
            'plane_deviation_span': (
                self._plane_deviation_span(deviations)
            ),
            'planarity_parameter': (
                self._planarity_parameter(deviations)
            ),
        }

    def get_results(
        self,
        mol,
        plane_atom_ids=None,
        deviation_atom_ids=None,
    ):
        """
        Calculate the planarity of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose planarity is to be calculated.

        plane_atom_ids : iterable of :class:`int`, optional
            The atom ids to use to define the plane of best fit.

        deviation_atom_ids : iterable of :class:`int`, optional
            The atom ids to use to calculate planarity.

        Returns
        -------
        :class:`.PlanarityResults`
            The planarity measures of the molecule.

        """

        return PlanarityResults(self.calculate(
            mol=mol,
            plane_atom_ids=plane_atom_ids,
            deviation_atom_ids=deviation_atom_ids,
        ))
