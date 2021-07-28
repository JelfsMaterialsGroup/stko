"""
Orca Extractor
=============

#. :class:`.OrcaExtractor`

Class to extract properties from Orca output.

"""

import re
from .extractor import Extractor


class OrcaExtractor(Extractor):
    """
    Extracts properties from Orca 4.2 output files.

    Limited to final single point energy for now.

    Attributes
    ----------
    output_file : :class:`str`
        Output file to extract properties from.

    output_lines : :class:`list` : :class:`str`
        :class:`list` of all lines in as :class:`str` in the output
        file.

    total_energy : :class:`float`
        The total energy in the :attr:`output_file` as
        :class:`float`. The energy is in units of a.u..

    """

    def _extract_values(self):
        """
        Extract all properties from Orca output file.

        Returns
        -------
        None : :class:`NoneType`

        """

        for i, line in enumerate(self.output_lines):
            if self._check_line(line, 'total_energy'):
                self._extract_total_energy(line)

    def _properties_dict(self):

        return {'total_energy': 'FINAL SINGLE POINT ENERGY'}

    def _extract_total_energy(self, line):
        """
        Updates :attr:`total_energy`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        Returns
        -------
        None : :class:`NoneType`

        """

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip()).group(0)
        self.total_energy = float(string)
