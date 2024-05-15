from collections import abc
from pathlib import Path

from stko._internal.calculators.extractors.orca_extractor import OrcaExtractor


class OrcaResults:
    """Results class containing molecule Orca properties."""

    def __init__(
        self,
        generator: abc.Iterable,
        output_file: Path | str,
        extractor: type = OrcaExtractor,
    ) -> None:
        # Run calculation.
        next(generator)  # type: ignore[call-overload]
        self._extractor = extractor(output_file=output_file)

    def get_total_energy(self) -> tuple[float, str]:
        return (self._extractor.total_energy, "a.u.")
