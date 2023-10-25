"""
XTB Extractor
=============

#. :class:`.XTBExtractor`

Class to extract properties from xTB output.

"""

import re

from .extractor import Extractor


class XTBExtractor(Extractor):
    """
    Extracts properties from xTB output files.

    All formatting based on the 190418 version of xTB.

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

    homo_lumo_gap : :class:`float`
        The HOMO-LUMO gap in the :attr:`output_file` as
        :class:`float`. The gap is in units of eV.

    fermi_level : :class:`float`
        The Fermi level in the :attr:`output_file` as
        :class:`float` in units of eV.

    qonly_dipole_moment : :class:`list`
        Components of the Q only dipole moment in units
        of Debye in :class:`list` of the form
        ``[x, y, z]``.

    full_dipole_moment : :class:`list`
        Components of the full dipole moment in units
        of Debye in :class:`list` of the form
        ``[x, y, z, total]``.

    qonly_quadrupole_moment : :class:`list`
        Components of the Q only traceless quadrupole moment in units
        of Debye in :class:`list` of the form
        ``[xx, xy, xy, xz, yz, zz]``.

    qdip_quadrupole_moment : :class:`list`
        Components of the Q+Dip traceless quadrupole moment in units of
        Debye in :class:`list` of the form
        ``[xx, xy, xy, xz, yz, zz]``.

    full_quadrupole_moment : :class:`list`
        Components of the full traceless quadrupole moment in units of
        Debye in :class:`list` of the form
        ``[xx, xy, xy, xz, yz, zz]``.

    homo_lumo_occ : :class:`dict`
        :class:`dict` of :class:`list` containing the orbital number,
        energy in eV and occupation of the HOMO and LUMO orbitals in
        the :attr:`output_file`.

    total_free_energy : :class:`float`
        The total free energy in the :attr:`output_file` as
        :class:`float`. The free energy is in units of a.u. and
        calculated at 298.15K.

    frequencies : :class:`list`
        :class:`list` of the vibrational frequencies as :class:`float`
        in the :attr:`output_file`. Vibrational frequencies are in
        units of wavenumber and calculated at 298.15K.

    ionisation_potential : :class:`float`
        The vertical ionisation potential in the :attr:`output_file` as
        :class:`float`. Corresponds to the delta SCC IP.

    electron_affinity : :class:`float`
        The vertical electron affinity in the :attr:`output_file` as
        :class:`float`. Corresponds to the delta SCC EA.

    """

    def _extract_values(self):

        for i, line in enumerate(self.output_lines):
            if self._check_line(line, 'total_energy'):
                self._extract_total_energy(line)
            elif self._check_line(line, 'homo_lumo_gap'):
                self._extract_homo_lumo_gap(line)
            elif self._check_line(line, 'fermi_level'):
                self._extract_fermi_level(line)
            elif self._check_line(line, 'dipole_moment'):
                self._extract_qonly_dipole_moment(i)
                self._extract_full_dipole_moment(i)
            elif self._check_line(line, 'quadrupole_moment'):
                self._extract_qonly_quadrupole_moment(i)
                self._extract_qdip_quadrupole_moment(i)
                self._extract_full_quadrupole_moment(i)
            elif self._check_line(line, 'homo_lumo_occ_HOMO'):
                self.homo_lumo_occ = {}
                self._extract_homo_lumo_occ(line, 'HOMO')
            elif self._check_line(line, 'homo_lumo_occ_LUMO'):
                self._extract_homo_lumo_occ(line, 'LUMO')
            elif self._check_line(line, 'total_free_energy'):
                self._extract_total_free_energy(line)
            elif self._check_line(line, 'ionisation_potential'):
                self._extract_ionisation_potential(line)
            elif self._check_line(line, 'electron_affinity'):
                self._extract_electron_affinity(line)

        # Frequency formatting requires loop through full file.
        self._extract_frequencies()

    def _properties_dict(self):

        return {
            'total_energy': '          | TOTAL ENERGY  ',
            'homo_lumo_gap': '          | HOMO-LUMO GAP   ',
            'fermi_level': '             Fermi-level        ',
            'dipole_moment': 'molecular dipole:',
            'quadrupole_moment': 'molecular quadrupole (traceless):',
            'homo_lumo_occ_HOMO': '(HOMO)',
            'homo_lumo_occ_LUMO': '(LUMO)',
            'total_free_energy': '          | TOTAL FREE ENERGY  ',
            'ionisation_potential': 'delta SCC IP (eV)',
            'electron_affinity': 'delta SCC EA (eV)',
        }

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

    def _extract_homo_lumo_gap(self, line):
        """
        Updates :attr:`homo_lumo_gap`.

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
        self.homo_lumo_gap = float(string)

    def _extract_fermi_level(self, line):
        """
        Updates :attr:`fermi_level`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        Returns
        -------
        None : :class:`NoneType`

        """

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        part2 = line.split('Eh')
        string = nums.search(part2[1].rstrip()).group(0)
        self.fermi_level = float(string)

    def _extract_qonly_dipole_moment(self, index):
        """
        Updates :attr:`qonly_dipole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+2].rstrip()

        if 'q only:' in sample_set:
            self.qonly_dipole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_full_dipole_moment(self, index):
        """
        Updates :attr:`full_dipole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+3].rstrip()

        if 'full:' in sample_set:
            self.full_dipole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_qonly_quadrupole_moment(self, index):
        """
        Updates :attr:`qonly_quadrupole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+2].rstrip()

        if 'q only:' in sample_set:
            self.qonly_quadrupole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_qdip_quadrupole_moment(self, index):
        """
        Updates :attr:`qdip_quadrupole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+3].rstrip()

        if 'q+dip:' in sample_set:
            self.qdip_quadrupole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_full_quadrupole_moment(self, index):
        """
        Updates :attr:`full_quadrupole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+4].rstrip()

        if 'full:' in sample_set:
            self.full_quadrupole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_homo_lumo_occ(self, line, orbital):
        """
        Updates :attr:`homo_lumo_occ`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        orbital : :class:`str`
            Can be 'HOMO' or 'LUMO'.

        Returns
        -------
        None : :class:`NoneType`

        """

        if orbital == 'HOMO':
            split_line = [i for i in line.rstrip().split(' ') if i]
            # The line is:
            #   Number, occupation, energy (Ha), energy (ev), label
            # Extract:
            #   Number, occupation, energy (eV)
            orbital_val = [
                int(split_line[0]),
                float(split_line[1]),
                float(split_line[3])
            ]
        elif orbital == 'LUMO':
            split_line = [i for i in line.rstrip().split(' ') if i]
            # The line is:
            #   Number, energy (Ha), energy (ev), label
            # Extract:
            #   Number, occupation (zero), energy (eV)
            orbital_val = [
                int(split_line[0]),
                0,
                float(split_line[2])
            ]

        self.homo_lumo_occ[orbital] = orbital_val

    def _extract_total_free_energy(self, line):
        """
        Updates :attr:`total_free_energy`.

        Returns
        -------
        None : :class:`NoneType`

        """

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip()).group(0)
        self.total_free_energy = float(string)

    def _extract_frequencies(self):
        """
        Updates :attr:`frequencies`.

        Returns
        -------
        None : :class:`NoneType`

        """

        test = '|               Frequency Printout                |'

        # Use a switch to make sure we are extracting values after the
        # final property readout.
        switch = False

        frequencies = []
        for i, line in enumerate(self.output_lines):
            if test in line:
                # Turn on reading as final frequency printout has
                # begun.
                switch = True
            if ' reduced masses (amu)' in line:
                # Turn off reading as frequency section is done.
                switch = False
            if 'eigval :' in line and switch is True:
                samp = line.rstrip().split(':')[1].split(' ')
                split_line = [i for i in samp if i]
                for freq in split_line:
                    frequencies.append(freq)

        self.frequencies = [float(i) for i in frequencies]

    def _extract_ionisation_potential(self, line):
        """
        Updates :attr:`ionisation_potential`.

        Returns
        -------
        None : :class:`NoneType`

        """

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip()).group(0)
        self.ionisation_potential = float(string)

    def _extract_electron_affinity(self, line):
        """
        Updates :attr:`electron_affinity`.

        Returns
        -------
        None : :class:`NoneType`

        """

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip()).group(0)
        self.electron_affinity = float(string)
