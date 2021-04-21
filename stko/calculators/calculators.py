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

    def calculate(self, mol):
        """
        Perform calculation on `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` to calculate properties of.

        Yields
        ------
        :class:`function`
            The function to perform the calculation.

        """

        raise NotImplementedError()

    def get_results(self, mol):
        """
        Get results of calculation on `mol`.

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
