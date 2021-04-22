"""
ShapeResults
=============

"""

from .results import Results


class ShapeResults(Results):
    """
    Results class containing molecule energy.

    """

    def __init__(self, generator):

        self._values = next(generator)
        print(self._values)
