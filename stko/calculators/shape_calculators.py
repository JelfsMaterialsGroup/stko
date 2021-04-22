"""
Shape Calculators
=================

#. :class:`.ShapeCalculator`

"""

import logging
from rdkit.Chem import AllChem as rdkit

from .calculators import Calculator
from .results import ShapeResults


logger = logging.getLogger(__name__)


class ShapeCalculator(Calculator):
    """
    Calculates shape measures of a molecule.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        FULLLL

    """

    def calculate(self, mol):
        yield

    def get_results(self, mol):
        """
        Calculate the shape measures of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`.ShapeResults`
            The shape measures of the molecule.

        """

        return ShapeResults(
            generator=self.calculate(mol),
        )
