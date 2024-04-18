from collections import abc


class PlanarityResults:
    """Results class containing molecule planarity measures."""

    def __init__(self, generator: abc.Iterable) -> None:
        self._values = next(generator)  # type: ignore[call-overload]

    def get_planarity_parameter(self) -> float:
        return self._values["planarity_parameter"]

    def get_plane_deviation(self) -> float:
        return self._values["plane_deviation"]

    def get_plane_deviation_span(self) -> float:
        return self._values["plane_deviation_span"]
