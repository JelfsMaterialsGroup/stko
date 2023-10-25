"""
Shape Results
=============

#. :class:`.ShapeResults`

Results class for extracting molecular shape properties.

"""

from .results import Results


class ShapeResults(Results):
    """
    Results class containing molecule shape measures.

    """

    def __init__(self, generator):

        self._values = next(generator)

    def get_pmi1(self):
        return self._values['pmi1']

    def get_pmi2(self):
        return self._values['pmi2']

    def get_pmi3(self):
        return self._values['pmi3']

    def get_npr1(self):
        return self._values['npr1']

    def get_npr2(self):
        return self._values['npr2']

    def get_asphericity(self):
        return self._values['asphericity']

    def get_eccentricity(self):
        return self._values['eccentricity']

    def get_inertial_shape_factor(self):
        return self._values['inertialshapefactor']

    def get_radius_of_gyration(self):
        return self._values['radiusofgyration']

    def get_spherocity_index(self):
        return self._values['spherocityindex']
