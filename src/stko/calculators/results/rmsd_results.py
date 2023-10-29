"""
RMSD Results
=============

Results class for extracting RMSD of two molecules.

"""


class RmsdResults:
    """
    Results class containing RMSD measures.

    """

    def __init__(self, generator):
        self._value = next(generator)

    def get_rmsd(self):
        return self._value
