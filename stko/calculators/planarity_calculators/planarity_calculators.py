"""
Planarity Calculators
=====================

#. :class:`.PlanarityCalculator`

Methods to calculate planarity measures of a molecule.

"""

import logging

from ..calculators import Calculator
from ..results import PlanarityResults

logger = logging.getLogger(__name__)


class PlanarityCalculator(Calculator):
    """
    Calculates measures of planarity of a molecule.

    Measures based on plane deviation from Angew. paper [1]_ and a
    ChemRxiv paper [2]_.

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
        ADDD

    References
    ----------
    .. [1] https://onlinelibrary.wiley.com/doi/10.1002/anie.202106721

    .. [2] https://chemrxiv.org/engage/chemrxiv/article-details/
    60c73cbf9abda2e0c5f8b5c6

    """

    def _get_plane_of_best_fit(self, mol, plane_atom_ids):
        raise NotImplementedError()

    def _plane_deviation(
        self,
        mol,
        plane_atom_ids=None,
        deviation_atom_ids=None,
    ):
        raise NotImplementedError()

    def _plane_deviation_span(
        self,
        mol,
        plane_atom_ids=None,
        deviation_atom_ids=None,
    ):
        raise NotImplementedError()

    def _planarity_parameter(
        self,
        mol,
        plane_atom_ids=None,
        deviation_atom_ids=None,
    ):
        raise NotImplementedError()

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

        yield self._plane_deviation(
            mol=mol,
            plane_atom_ids=plane_atom_ids,
            deviation_atom_ids=deviation_atom_ids,
        )
        yield self._plane_deviation_span(
            mol=mol,
            plane_atom_ids=plane_atom_ids,
            deviation_atom_ids=deviation_atom_ids,
        )
        yield self._planarity_parameter(
            mol=mol,
            plane_atom_ids=plane_atom_ids,
            deviation_atom_ids=deviation_atom_ids,
        )

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
