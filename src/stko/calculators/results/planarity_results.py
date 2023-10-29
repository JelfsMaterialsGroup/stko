"""
Planarity Results
=================

Results class for extracting molecular planarity measures.

"""


class PlanarityResults:
    """
    Results class containing molecule planarity measures.

    """

    def __init__(self, generator):
        self._values = next(generator)

    def get_planarity_parameter(self):
        return self._values["planarity_parameter"]

    def get_plane_deviation(self):
        return self._values["plane_deviation"]

    def get_plane_deviation_span(self):
        return self._values["plane_deviation_span"]
