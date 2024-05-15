import logging
from collections import abc

import numpy as np
import stk

from stko._internal.calculators.results.planarity_results import (
    PlanarityResults,
)

logger = logging.getLogger(__name__)


class PlanarityCalculator:
    """Calculates measures of planarity of a molecule.

    Measures based on plane deviation from Angew. paper [1]_ and a
    ChemRxiv paper [2]_.

    Plane deviation: sum of the shortest distance to the plane of best
    fit of all deviation atoms (sum abs(d_i)).

    Plane deviation span: d_max - d_min (SDP in [2]_)

    Planarity parameter: defined as
    sqrt((1/num_atoms) * (sum d_i ** 2)) (MPP in [2]_)

    Examples:
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

    References:
        .. [1] https://onlinelibrary.wiley.com/doi/10.1002/anie.202106721

        .. [2] https://link.springer.com/article/10.1007/s00894-021-04884-0

    """

    def _get_plane_of_best_fit(
        self,
        mol: stk.Molecule,
        plane_atom_ids: abc.Iterable[int],
    ) -> np.ndarray:
        centroid = mol.get_centroid(atom_ids=plane_atom_ids)
        normal = mol.get_plane_normal(atom_ids=plane_atom_ids)
        # Plane of equation ax + by + cz = d.
        return np.append(normal, np.sum(normal * centroid))

    def _shortest_distance_to_plane(
        self,
        plane: np.ndarray,
        point: np.ndarray,
    ) -> float:
        """Calculate the perpendicular distance from a point and a plane."""
        top = (
            plane[0] * point[0]
            + plane[1] * point[1]
            + plane[2] * point[2]
            - plane[3]
        )
        bottom = np.sqrt(plane[0] ** 2 + plane[1] ** 2 + plane[2] ** 2)
        return top / bottom

    def _calculate_deviations(
        self,
        mol: stk.Molecule,
        atom_plane: np.ndarray,
        deviation_atom_ids: abc.Iterable[int],
    ) -> list[float]:
        return [
            self._shortest_distance_to_plane(
                plane=atom_plane,
                point=next(
                    mol.get_atomic_positions(atom_ids=i.get_id()),
                ),
            )
            for i in mol.get_atoms()
            if i.get_id() in deviation_atom_ids
        ]

    def _plane_deviation(self, deviations: list[float]) -> float:
        deviations = [abs(i) for i in deviations]
        return sum(deviations)

    def _plane_deviation_span(self, deviations: list[float]) -> float:
        return max(deviations) - min(deviations)

    def _planarity_parameter(self, deviations: list[float]) -> float:
        deviations = [abs(i) for i in deviations]
        num_atoms = len(deviations)
        inv_num_atoms = 1 / num_atoms
        sum_squared = sum(i**2 for i in deviations)
        return np.sqrt(inv_num_atoms * sum_squared)

    def calculate(
        self,
        mol: stk.Molecule,
        plane_atom_ids: abc.Iterable[int] | None = None,
        deviation_atom_ids: abc.Iterable[int] | None = None,
    ) -> abc.Iterable[dict]:
        """Perform calculation on `mol`.

        Parameters:
            mol:
                The :class:`stk.Molecule` whose planarity is to be calculated.

            plane_atom_ids:
                The atom ids to use to define the plane of best fit.

            deviation_atom_ids:
                The atom ids to use to calculate planarity.

        Yields:
            Dictionary of results.

        """
        if plane_atom_ids is None:
            plane_atom_ids = list(range(len(list(mol.get_atoms()))))
        else:
            plane_atom_ids = list(plane_atom_ids)

        if deviation_atom_ids is None:
            deviation_atom_ids = list(range(len(list(mol.get_atoms()))))
        else:
            deviation_atom_ids = list(deviation_atom_ids)

        atom_plane = self._get_plane_of_best_fit(mol, plane_atom_ids)
        deviations = self._calculate_deviations(
            mol=mol,
            atom_plane=atom_plane,
            deviation_atom_ids=deviation_atom_ids,
        )
        yield {
            "plane_deviation": (self._plane_deviation(deviations)),
            "plane_deviation_span": (self._plane_deviation_span(deviations)),
            "planarity_parameter": (self._planarity_parameter(deviations)),
        }

    def get_results(
        self,
        mol: stk.Molecule,
        plane_atom_ids: abc.Iterable[int] | None = None,
        deviation_atom_ids: abc.Iterable[int] | None = None,
    ) -> PlanarityResults:
        """Calculate the planarity of `mol`.

        Parameters:
            mol:
                The :class:`stk.Molecule` whose planarity is to be calculated.

            plane_atom_ids:
                The atom ids to use to define the plane of best fit.

            deviation_atom_ids:
                The atom ids to use to calculate planarity.

        Returns:
            The planarity measures of the molecule.

        """
        return PlanarityResults(
            self.calculate(
                mol=mol,
                plane_atom_ids=plane_atom_ids,
                deviation_atom_ids=deviation_atom_ids,
            )
        )
