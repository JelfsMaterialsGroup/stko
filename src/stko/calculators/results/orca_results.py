"""
Orca Results
============

#. :class:`.OrcaResults`

Results class for the output of Orca.

"""

from .results import Results
from ..extractors import OrcaExtractor


class OrcaResults(Results):
    """
    Results class containing molecule Orca properties.

    """

    def __init__(
        self,
        generator,
        output_file,
        extractor=OrcaExtractor,
    ):

        # Run calculation.
        next(generator)
        self._extractor = extractor(output_file=output_file)

    def get_total_energy(self):
        return (self._extractor.total_energy, 'a.u.')
