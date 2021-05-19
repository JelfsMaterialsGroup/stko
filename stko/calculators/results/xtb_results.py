"""
XTB Results
===========

#. :class:`.XTBResults`

Results class for the output of XTB.

"""

from .results import Results
from ...utilities import XTBExtractor


class XTBResults(Results):
    """
    Results class containing molecule xTB properties.

    """

    def __init__(self, generator, output_file, extractor=XTBExtractor):

        # Run calculation.
        next(generator)
        self._extractor = extractor(output_file=output_file)

    def get_total_energy(self):
        return (self._extractor.total_energy, 'a.u.')

    def get_homo_lumo_gap(self):
        return (self._extractor.homo_lumo_gap, 'eV')

    def get_fermi_level(self):
        return (self._extractor.fermi_level, 'eV')

    def get_homo_lumo_orbitals(self):
        return (self._extractor.homo_lumo_occ, 'eV')

    def get_qonly_dipole_moments(self):
        return (self._extractor.qonly_dipole_moment, 'Debye')

    def get_full_dipole_moments(self):
        return (self._extractor.full_dipole_moment, 'Debye')

    def get_qonly_quadrupole_moments(self):
        return (self._extractor.qonly_quadrupole_moment, 'Debye')

    def get_qdip_quadrupole_moments(self):
        return (self._extractor.qdip_quadrupole_moment, 'Debye')

    def get_full_quadrupole_moments(self):
        return (self._extractor.full_quadrupole_moment, 'Debye')

    def get_total_free_energy(self):
        try:
            return (self._extractor.total_free_energy, 'a.u.')
        except AttributeError:
            raise AttributeError(
                'Frequency, hessian and thermo calculations not '
                'performed to extract this property.'
            )

    def get_frequencies(self):
        try:
            return (self._extractor.frequencies, 'wavenumber')
        except AttributeError:
            raise AttributeError(
                'Frequency, hessian and thermo calculations not '
                'performed to extract this property.'
            )
