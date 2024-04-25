from collections import abc


class ShapeResults:
    """Results class containing molecule shape measures."""

    def __init__(self, generator: abc.Iterable) -> None:
        self._values = next(generator)  # type: ignore[call-overload]

    def get_pmi1(self) -> float:
        return self._values["pmi1"]

    def get_pmi2(self) -> float:
        return self._values["pmi2"]

    def get_pmi3(self) -> float:
        return self._values["pmi3"]

    def get_npr1(self) -> float:
        return self._values["npr1"]

    def get_npr2(self) -> float:
        return self._values["npr2"]

    def get_asphericity(self) -> float:
        return self._values["asphericity"]

    def get_eccentricity(self) -> float:
        return self._values["eccentricity"]

    def get_inertial_shape_factor(self) -> float:
        return self._values["inertialshapefactor"]

    def get_radius_of_gyration(self) -> float:
        return self._values["radiusofgyration"]

    def get_spherocity_index(self) -> float:
        return self._values["spherocityindex"]
