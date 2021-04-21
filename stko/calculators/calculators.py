"""
Calculator
==========

"""

from abc import ABC


class Calculator(ABC):
    """
    An abstract base class for calculators.
    ...

    """

    def get_results(self, mol):
        """
        Perform calculation on `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` to calculate properties of.

        Returns
        -------
        :class:`.Results`
            The result of the calculation.

        """

        raise NotImplementedError()
