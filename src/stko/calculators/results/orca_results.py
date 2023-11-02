import typing

from stko.calculators.extractors.orca_extractor import OrcaExtractor


class OrcaResults:
    """
    Results class containing molecule Orca properties.

    """

    def __init__(
        self,
        generator: typing.Iterable,
        output_file: str,
        extractor: type = OrcaExtractor,
    ):
        # Run calculation.
        next(generator)  # type: ignore[call-overload]
        self._extractor = extractor(output_file=output_file)

    def get_total_energy(self) -> tuple[float, str]:
        return (self._extractor.total_energy, "a.u.")
