"""
TorsionResults
==============

"""

from .results import Results


class TorsionResults(Results):
    """
    Results class containing molecule torsions.

    """

    def __init__(self, generator):
        raise NotImplementedError()


class ConstructedMoleculeTorsionResults(TorsionResults):
    """
    Results class containing molecule torsions.

    """

    def __init__(self, generator):
        raise NotImplementedError()
