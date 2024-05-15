from dataclasses import dataclass

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    metal_atom_distances: dict[tuple[int, int], float]
    metal_centroid_angles: dict[tuple[int, int], float]
    min_centoid_distance: float
    avg_centoid_distance: tuple[float, float]
    radius_gyration: float
    min_atom_atom_distance: float
    max_diameter: float
    bonds: dict[tuple[str, str], list[float]]
    angles: dict[tuple[str, str, str], list[float]]
    torsions: dict[tuple[str, str, str, str], list[float]]
    name: str
