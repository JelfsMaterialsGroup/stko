"""
EnergyResults
=============

"""

from .results import Results


class EnergyResults(Results):
    """
    Results class containing molecule energy.

    """

    def __init__(self, value, unit_string):

        self._value = value
        self._unit_string = unit_string

    def get_energy(self):
        return self._value

    def get_unit_string(self):
        return self._unit_string
