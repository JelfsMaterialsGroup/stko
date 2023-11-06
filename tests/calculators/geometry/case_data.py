import stk


class CaseData:
    def __init__(
        self,
        molecule: stk.Molecule,
        metal_atom_distances: dict[tuple[int, int], float],
        metal_centroid_angles: dict[tuple[int, int], float],
        min_centoid_distance: float,
        avg_centoid_distance: tuple[float, float],
        radius_gyration: float,
        min_atom_atom_distance: float,
        max_diameter: float,
        bonds: dict[tuple[str, str], list[float]],
        angles: dict[tuple[str, str, str], list[float]],
        torsions: dict[tuple[str, str, str, str], list[float]],
        name: str,
    ) -> None:
        self.molecule = molecule
        self.metal_atom_distances = metal_atom_distances
        self.metal_centroid_angles = metal_centroid_angles
        self.min_centoid_distance = min_centoid_distance
        self.avg_centoid_distance = avg_centoid_distance
        self.radius_gyration = radius_gyration
        self.min_atom_atom_distance = min_atom_atom_distance
        self.max_diameter = max_diameter
        self.bonds = bonds
        self.angles = angles
        self.torsions = torsions
        self.name = name
