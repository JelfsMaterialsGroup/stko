import re
from pathlib import Path

from stko._internal.calculators.extractors.utilities import check_line


class XTBExtractor:
    """Extracts properties from xTB output files.

    All formatting based on the 190418 version of xTB.

    Parameters:
        output_file:
            Output file to extract properties from.

    Attributes:
        output_file:
            Output file to extract properties from.

        output_lines:
            List of all lines in as string in the output file.

        total_energy:
            The total energy in the :attr:`output_file`. The energy is
            in units of a.u..

        homo_lumo_gap:
            The HOMO-LUMO gap in the :attr:`output_file`. The gap is
            in units of eV.

        fermi_level:
            The Fermi level in the :attr:`output_file` in units of eV.

        qonly_dipole_moment:
            Components of the Q only dipole moment in units
            of Debye in List of the form ``[x, y, z]``.

        full_dipole_moment:
            Components of the full dipole moment in units
            of Debye in List of the form
            ``[x, y, z, total]``.

        qonly_quadrupole_moment:
            Components of the Q only traceless quadrupole moment in units
            of Debye in List of the form
            ``[xx, xy, xy, xz, yz, zz]``.

        qdip_quadrupole_moment:
            Components of the Q+Dip traceless quadrupole moment in units of
            Debye in List of the form
            ``[xx, xy, xy, xz, yz, zz]``.

        full_quadrupole_moment:
            Components of the full traceless quadrupole moment in units of
            Debye in List of the form
            ``[xx, xy, xy, xz, yz, zz]``.

        homo_lumo_occ:
            Dictionary of List containing the orbital number,
            energy in eV and occupation of the HOMO and LUMO orbitals in
            the :attr:`output_file`.

        total_free_energy:
            The total free energy in the :attr:`output_file`.
            The free energy is in units of a.u. and calculated at 298.15K.

        frequencies:
            List of the vibrational frequencies in the :attr:`output_file`.
            Vibrational frequencies are in units of wavenumber and
            calculated at 298.15K.

        ionisation_potential:
            The vertical ionisation potential in the :attr:`output_file`.
            Corresponds to the delta SCC IP.

        electron_affinity:
            The vertical electron affinity in the :attr:`output_file`.
            Corresponds to the delta SCC EA.

        total_sasa:
            The solvent-accessible surface area of the molecule from xtb.

    Examples:
        .. code-block:: python

            import stko

            data = stko.XTBExtractor(output_file)
            print(data.total_energy)
            print(data.homo_lumo_gap)

    """

    def __init__(self, output_file: Path | str) -> None:
        self.output_file = Path(output_file)
        # Explictly set encoding to UTF-8 because default encoding on
        # Windows will fail to read the file otherwise.
        with self.output_file.open(encoding="UTF-8") as f:
            self.output_lines = f.readlines()

        self._extract_values()

    def _extract_values(self) -> None:  # noqa: C901
        for i, line in enumerate(self.output_lines):
            if check_line(line, "total_energy", self._properties_dict()):
                self._extract_total_energy(line)
            elif check_line(line, "homo_lumo_gap", self._properties_dict()):
                self._extract_homo_lumo_gap(line)
            elif check_line(line, "fermi_level", self._properties_dict()):
                self._extract_fermi_level(line)
            elif check_line(line, "dipole_moment", self._properties_dict()):
                self._extract_qonly_dipole_moment(i)
                self._extract_full_dipole_moment(i)
            elif check_line(
                line, "quadrupole_moment", self._properties_dict()
            ):
                self._extract_qonly_quadrupole_moment(i)
                self._extract_qdip_quadrupole_moment(i)
                self._extract_full_quadrupole_moment(i)
            elif check_line(
                line, "homo_lumo_occ_HOMO", self._properties_dict()
            ):
                self.homo_lumo_occ: dict[str, list[float]] = {}
                self._extract_homo_lumo_occ(line, "HOMO")
            elif check_line(
                line, "homo_lumo_occ_LUMO", self._properties_dict()
            ):
                self._extract_homo_lumo_occ(line, "LUMO")
            elif check_line(
                line, "total_free_energy", self._properties_dict()
            ):
                self._extract_total_free_energy(line)
            elif check_line(
                line, "ionisation_potential", self._properties_dict()
            ):
                self._extract_ionisation_potential(line)
            elif check_line(
                line, "electron_affinity", self._properties_dict()
            ):
                self._extract_electron_affinity(line)
            elif check_line(line, "total_sasa", self._properties_dict()):
                self._extract_total_sasa(line)

        # Frequency formatting requires loop through full file.
        self._extract_frequencies()

    def _properties_dict(self) -> dict[str, str]:
        return {
            "total_energy": "          | TOTAL ENERGY  ",
            "homo_lumo_gap": "          | HOMO-LUMO GAP   ",
            "fermi_level": "             Fermi-level        ",
            "dipole_moment": "molecular dipole:",
            "quadrupole_moment": "molecular quadrupole (traceless):",
            "homo_lumo_occ_HOMO": "(HOMO)",
            "homo_lumo_occ_LUMO": "(LUMO)",
            "total_free_energy": "          | TOTAL FREE ENERGY  ",
            "ionisation_potential": "delta SCC IP (eV)",
            "electron_affinity": "delta SCC EA (eV)",
            "total_sasa": "total SASA /",
        }

    def _extract_total_energy(self, line: str) -> None:
        """Updates :attr:`total_energy`.

        Parameters:
            line:
                Line of output file to extract property from.

        """
        # Use regex to match to numbers.
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip())
        self.total_energy = float(string.group(0))  # type: ignore[union-attr]

    def _extract_homo_lumo_gap(self, line: str) -> None:
        """Updates :attr:`homo_lumo_gap`.

        Parameters:
            line:
                Line of output file to extract property from.

        """
        # Use regex to match to numbers.
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip())
        self.homo_lumo_gap = float(string.group(0))  # type: ignore[union-attr]

    def _extract_fermi_level(self, line: str) -> None:
        """Updates :attr:`fermi_level`.

        Parameters:
            line:
                Line of output file to extract property from.

        """
        # Use regex to match to numbers.
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        part2 = line.split("Eh")
        string = nums.search(part2[1].rstrip())
        self.fermi_level = float(string.group(0))  # type: ignore[union-attr]

    def _extract_qonly_dipole_moment(self, index: int) -> None:
        """Updates :attr:`qonly_dipole_moment`.

        Parameters:
            index:
                Index of line in :attr:`output_lines`.

        """
        sample_set = self.output_lines[index + 2].rstrip()

        if "q only:" in sample_set:
            self.qonly_dipole_moment = [
                float(i) for i in sample_set.split(":")[1].split(" ") if i
            ]

    def _extract_full_dipole_moment(self, index: int) -> None:
        """Updates :attr:`full_dipole_moment`.

        Parameters:
            index:
                Index of line in :attr:`output_lines`.

        """
        sample_set = self.output_lines[index + 3].rstrip()

        if "full:" in sample_set:
            self.full_dipole_moment = [
                float(i) for i in sample_set.split(":")[1].split(" ") if i
            ]

    def _extract_qonly_quadrupole_moment(self, index: int) -> None:
        """Updates :attr:`qonly_quadrupole_moment`.

        Parameters:
            index:
                Index of line in :attr:`output_lines`.

        """
        sample_set = self.output_lines[index + 2].rstrip()

        if "q only:" in sample_set:
            self.qonly_quadrupole_moment = [
                float(i) for i in sample_set.split(":")[1].split(" ") if i
            ]

    def _extract_qdip_quadrupole_moment(self, index: int) -> None:
        """Updates :attr:`qdip_quadrupole_moment`.

        Parameters:
            index:
                Index of line in :attr:`output_lines`.

        """
        sample_set = self.output_lines[index + 3].rstrip()

        if "q+dip:" in sample_set:
            self.qdip_quadrupole_moment = [
                float(i) for i in sample_set.split(":")[1].split(" ") if i
            ]

    def _extract_full_quadrupole_moment(self, index: int) -> None:
        """Updates :attr:`full_quadrupole_moment`.

        Parameters:
            index:
                Index of line in :attr:`output_lines`.

        """
        sample_set = self.output_lines[index + 4].rstrip()

        if "full:" in sample_set:
            self.full_quadrupole_moment = [
                float(i) for i in sample_set.split(":")[1].split(" ") if i
            ]

    def _extract_homo_lumo_occ(self, line: str, orbital: str) -> None:
        """Updates :attr:`homo_lumo_occ`.

        Parameters:
            line:
                Line of output file to extract property from.

            orbital:
                Can be 'HOMO' or 'LUMO'.

        """
        if orbital == "HOMO":
            split_line = [i for i in line.rstrip().split(" ") if i]
            # The line is: Number, occupation, energy (Ha), energy (ev), label
            # Extract: Number, occupation, energy (eV)
            orbital_val = [
                int(split_line[0]),
                float(split_line[1]),
                float(split_line[3]),
            ]
        elif orbital == "LUMO":
            split_line = [i for i in line.rstrip().split(" ") if i]
            # The line is: Number, energy (Ha), energy (ev), label
            # Extract: Number, occupation (zero), energy (eV)
            orbital_val = [int(split_line[0]), 0, float(split_line[2])]

        self.homo_lumo_occ[orbital] = orbital_val

    def _extract_total_free_energy(self, line: str) -> None:
        """Updates :attr:`total_free_energy`.

        Parameters:
            line:
                Line of output file to extract property from.

        """
        # Use regex to match to numbers.
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip())
        self.total_free_energy = float(
            string.group(0)  # type: ignore[union-attr]
        )

    def _extract_frequencies(self) -> None:
        """Updates :attr:`frequencies`."""
        test = "|               Frequency Printout                |"

        # Use a switch to make sure we are extracting values after the
        # final property readout.
        switch = False

        frequencies = []
        for _, line in enumerate(self.output_lines):
            if test in line:
                # Turn on reading as final frequency printout has
                # begun.
                switch = True
            if " reduced masses (amu)" in line:
                # Turn off reading as frequency section is done.
                switch = False
            if "eigval :" in line and switch is True:
                samp = line.rstrip().split(":")[1].split(" ")
                split_line = [i for i in samp if i]
                frequencies.extend(split_line)

        self.frequencies = [float(i) for i in frequencies]

    def _extract_ionisation_potential(self, line: str) -> None:
        """Updates :attr:`ionisation_potential`.

        Parameters:
            line:
                Line of output file to extract property from.

        """
        # Use regex to match to numbers.
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip())
        self.ionisation_potential = float(
            string.group(0)  # type: ignore[union-attr]
        )

    def _extract_electron_affinity(self, line: str) -> None:
        """Updates :attr:`electron_affinity`.

        Parameters:
            line:
                Line of output file to extract property from.

        """
        # Use regex to match to numbers.
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip())
        self.electron_affinity = float(
            string.group(0)  # type: ignore[union-attr]
        )

    def _extract_total_sasa(self, line: str) -> None:
        """Updates :attr:`total_sasa`.

        Parameters:
            line:
                Line of output file to extract property from.

        """
        # Use regex to match to numbers.
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip())
        self.total_sasa = float(string.group(0))  # type: ignore[union-attr]
