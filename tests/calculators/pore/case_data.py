import stk


class CaseData:
    def __init__(
        self,
        molecule: stk.Molecule,
        metal_atom_distances: dict[tuple[int, int], float],
        metal_centroid_angles: dict[tuple[int, int], float],
        min_centoid_distance: float,
        avg_centoid_distance: tuple[float, float],
        name: str,
    ) -> None:
        self.molecule = molecule
        self.metal_atom_distances = metal_atom_distances
        self.metal_centroid_angles = metal_centroid_angles
        self.min_centoid_distance = min_centoid_distance
        self.avg_centoid_distance = avg_centoid_distance
        self.name = name
