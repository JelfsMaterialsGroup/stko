"""
Planarity Results
=================

#. :class:`.PlanarityResults`

Results class for extracting molecular planarity measures.

"""

from .results import Results


class PlanarityResults(Results):
    """
    Results class containing molecule planarity measures.

    """

    def __init__(self, generator):

        self._values = next(generator)

    def get_planarity(self):
        raise NotImplementedError()
