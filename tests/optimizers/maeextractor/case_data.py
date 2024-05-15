from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class CaseData:
    run_name: str
    mae_path: str
    n: int
    energies: list[tuple[float, int]]
    min_energy: float
    path: str
    name: str
