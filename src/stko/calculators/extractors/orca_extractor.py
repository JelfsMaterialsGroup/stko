"""
Orca Extractor
==============

Class to extract properties from Orca output.

"""

import re

from stko.calculators.extractors.utilities import check_line


class OrcaExtractor:
    """
    Extracts properties from Orca 4.2 output files.

    Limited to final single point energy for now.

    Attributes:
        output_file:
            Output file to extract properties from.

        output_lines:
            :class:`list` of all lines in as :class:`str` in the output
            file.

        total_energy:
            The total energy in the :attr:`output_file` as
            :class:`float`. The energy is in units of a.u..

    Examples:

        .. code-block:: python

            import stko

            data = stko.OrcaExtractor(output_file)
            print(data.total_energy)


    """

    def __init__(self, output_file: str) -> None:
        """
        Initializes :class:`OrcaExtractor`

        Parameters:
            output_file:
                Output file to extract properties from.

        """

        self.output_file = output_file
        # Explictly set encoding to UTF-8 because default encoding on
        # Windows will fail to read the file otherwise.
        with open(self.output_file, "r", encoding="UTF-8") as f:
            self.output_lines = f.readlines()

        self._extract_values()

    def _extract_values(self):
        """
        Updates all properties by extracting from Orca output file.

        """

        for i, line in enumerate(self.output_lines):
            if check_line(line, "total_energy"):
                self._extract_total_energy(line)

    def _properties_dict(self) -> dict[str, str]:
        return {"total_energy": "FINAL SINGLE POINT ENERGY"}

    def _extract_total_energy(self, line: str):
        """
        Updates :attr:`total_energy`.

        Parameters:

            line:
            Line of output file to extract property from.

        """

        nms = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nms.search(line.rstrip()).group(0)  # type: ignore[union-attr]
        self.total_energy = float(string)
