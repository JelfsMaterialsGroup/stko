from collections import abc


class EnergyResults:
    """Results class containing molecule energy."""

    def __init__(self, generator: abc.Iterable, unit_string: str) -> None:
        self._value = next(generator)  # type: ignore[call-overload]
        self._unit_string = unit_string

    def get_energy(self) -> float:
        return self._value

    def get_unit_string(self) -> str:
        return self._unit_string
