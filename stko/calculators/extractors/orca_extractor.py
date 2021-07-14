"""
Orca Extractor
=============

#. :class:`.OrcaExtractor`

Class to extract properties from Orca output.

"""

import re


class OrcaExtractor:
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

    def __init__(self, output_file):
        """
        Initializes :class:`OrcaExtractor`

        Parameters
        ----------
        output_file : :class:`str`
            Output file to extract properties from.

        """

        self.output_file = output_file
        # Explictly set encoding to UTF-8 because default encoding on
        # Windows will fail to read the file otherwise.
        with open(self.output_file, 'r', encoding='UTF-8') as f:
            self.output_lines = f.readlines()

        self._extract_values()

    def _extract_values(self):
        """
        Extract all properties from xTB output file.

        Returns
        -------
        None : :class:`NoneType`

        """

        for i, line in enumerate(self.output_lines):
            if self._check_line(line, 'total_energy'):
                self._extract_total_energy(line)

    def _check_line(self, line, option):
        """
        Checks a line for a string based on option.

        All formatting based on the 190418 version of xTB.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to check.

        option : :class:`str`
            Define which property and string being checked for.
            Can be one of ``'total_energy'``.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if the desired string is present.

        """

        options = {
            'total_energy': 'FINAL SINGLE POINT ENERGY',
        }

        if options[option] in line:
            return True

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
