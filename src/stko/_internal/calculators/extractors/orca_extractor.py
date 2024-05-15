import re
from pathlib import Path

from stko._internal.calculators.extractors.utilities import check_line


class OrcaExtractor:
    """Extracts properties from Orca 4.2 output files.

    Limited to final single point energy for now.

    Parameters:
        output_file:
            Output file to extract properties from.


    Attributes:
        output_file:
            Output file to extract properties from.

        output_lines:
            List of all lines in as strings in the output file.

        total_energy:
            The total energy in the :attr:`output_file`.
            The energy is in units of a.u..

    Examples:
        .. code-block:: python

            import stko

            data = stko.OrcaExtractor(output_file)
            print(data.total_energy)


    """

    def __init__(self, output_file: Path | str) -> None:
        self.output_file = Path(output_file)
        # Explictly set encoding to UTF-8 because default encoding on
        # Windows will fail to read the file otherwise.
        with self.output_file.open(encoding="UTF-8") as f:
            self.output_lines = f.readlines()

        self._extract_values()

    def _extract_values(self) -> None:
        """Updates all properties by extracting from Orca output file."""
        for _, line in enumerate(self.output_lines):
            if check_line(line, "total_energy", self._properties_dict()):
                self._extract_total_energy(line)

    def _properties_dict(self) -> dict[str, str]:
        return {"total_energy": "FINAL SINGLE POINT ENERGY"}

    def _extract_total_energy(self, line: str) -> None:
        """Updates :attr:`total_energy`.

        Parameters:
            line: Line of output file to extract property from.

        """
        nms = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nms.search(line.rstrip()).group(0)  # type: ignore[union-attr]
        self.total_energy = float(string)
