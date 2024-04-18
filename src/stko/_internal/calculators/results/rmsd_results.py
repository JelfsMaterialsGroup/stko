class RmsdResults:
    """Results class containing RMSD measures."""

    def __init__(self, value: float) -> None:
        self._value = value

    def get_rmsd(self) -> float:
        return self._value
