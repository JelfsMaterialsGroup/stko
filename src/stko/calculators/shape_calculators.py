"""
Shape Calculators
=================

#. :class:`.ShapeCalculator`

Wrappers for shape calculators within the :mod:`rdkit` code.

"""

import logging
from rdkit.Chem import Descriptors3D as D3D

from .calculators import Calculator
from .results import ShapeResults


logger = logging.getLogger(__name__)


class ShapeCalculator(Calculator):
    """
    Calculates shape measures of a molecule.

    Uses :mod:`rdkit` 3D Descriptors [1]_ module to calculate all
    measures.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        bb1 = stk.BuildingBlock('C1CCCCC1')
        shape_calc = stko.ShapeCalculator()
        shape_results = shape_calc.get_results(bb1)
        eccentricity  = shape_results.get_eccentricity()

    References
    ----------
    .. [1]
    https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors3D.html

    """

    def calculate(self, mol):
        rdkit_mol = mol.to_rdkit_mol()
        results_dict = {
            'pmi1': D3D.PMI1(rdkit_mol),
            'pmi2': D3D.PMI2(rdkit_mol),
            'pmi3': D3D.PMI3(rdkit_mol),
            'npr1': D3D.NPR1(rdkit_mol),
            'npr2': D3D.NPR2(rdkit_mol),
            'asphericity': D3D.Asphericity(rdkit_mol),
            'eccentricity': D3D.Eccentricity(rdkit_mol),
            'inertialshapefactor': D3D.InertialShapeFactor(rdkit_mol),
            'radiusofgyration': D3D.RadiusOfGyration(rdkit_mol),
            'spherocityindex': D3D.SpherocityIndex(rdkit_mol),
        }

        yield results_dict

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

        return ShapeResults(self.calculate(mol))
