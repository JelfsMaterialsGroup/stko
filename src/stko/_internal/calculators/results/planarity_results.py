from collections import abc

import numpy as np


class PlanarityResults:
    """Results class containing molecule planarity measures."""

    def __init__(self, generator: abc.Iterable) -> None:
        self._values = next(generator)  # type: ignore[call-overload]

    def get_atom_plane(self) -> np.ndarray:
        return self._values["atom_plane"]

    def get_planarity_parameter(self) -> float:
        return float(self._values["planarity_parameter"])

    def get_plane_deviation(self) -> float:
        return float(self._values["plane_deviation"])

    def get_plane_deviation_span(self) -> float:
        return float(self._values["plane_deviation_span"])
