class CaseData:
    """A test case."""

    def __init__(
        self,
        run_name: str,
        mae_path: str,
        n: int,
        energies: list,
        min_energy: float,
        path: str,
        name: str,
    ) -> None:
        self.run_name = run_name
        self.n = n
        self.mae_path = mae_path
        self.energies = energies
        self.min_energy = min_energy
        self.path = path
        self.name = name
