"""
RMSD Results
=============

#. :class:`.RmsdResults`

Results class for extracting RMSD of two molecules.

"""

from stko.calculators.results.results import Results


class RmsdResults(Results):
    """
    Results class containing RMSD measures.

    """

    def __init__(self, generator):
        self._value = next(generator)

    def get_rmsd(self):
        return self._value
