from collections import abc
from pathlib import Path

from stko._internal.calculators.extractors.xtb_extractor import XTBExtractor


class XTBResults:
    """Results class containing molecule xTB properties."""

    def __init__(
        self,
        generator: abc.Generator,
        output_file: Path | str,
        extractor: type = XTBExtractor,
    ) -> None:
        # Run calculation.
        next(generator)  # type: ignore[call-overload]
        self._extractor = extractor(output_file=output_file)

    def get_total_energy(self) -> tuple[float, str]:
        return (self._extractor.total_energy, "a.u.")

    def get_homo_lumo_gap(self) -> tuple[float, str]:
        return (self._extractor.homo_lumo_gap, "eV")

    def get_fermi_level(self) -> tuple[float, str]:
        return (self._extractor.fermi_level, "eV")

    def get_homo_lumo_orbitals(self) -> tuple[float, str]:
        return (self._extractor.homo_lumo_occ, "eV")

    def get_qonly_dipole_moments(self) -> tuple[float, str]:
        return (self._extractor.qonly_dipole_moment, "Debye")

    def get_full_dipole_moments(self) -> tuple[float, str]:
        return (self._extractor.full_dipole_moment, "Debye")

    def get_qonly_quadrupole_moments(self) -> tuple[float, str]:
        return (self._extractor.qonly_quadrupole_moment, "Debye")

    def get_qdip_quadrupole_moments(self) -> tuple[float, str]:
        return (self._extractor.qdip_quadrupole_moment, "Debye")

    def get_full_quadrupole_moments(self) -> tuple[float, str]:
        return (self._extractor.full_quadrupole_moment, "Debye")

    def get_total_free_energy(self) -> tuple[float, str]:
        try:
            return (self._extractor.total_free_energy, "a.u.")
        except AttributeError as exc:
            msg = (
                "Frequency, hessian and thermo calculations not "
                "performed to extract this property."
            )
            raise AttributeError(msg) from exc

    def get_frequencies(self) -> tuple[float, str]:
        try:
            return (self._extractor.frequencies, "wavenumber")
        except AttributeError as exc:
            msg = (
                "Frequency, hessian and thermo calculations not "
                "performed to extract this property."
            )
            raise AttributeError(msg) from exc

    def get_ionisation_potential(self) -> tuple[float, str]:
        return (self._extractor.ionisation_potential, "eV")

    def get_electron_affinity(self) -> tuple[float, str]:
        return (self._extractor.electron_affinity, "eV")

    def get_total_sasa(self) -> tuple[float, str]:
        return (self._extractor.total_sasa, "A^2")
