import stk


class CaseData:
    """
    A test case.

    Attributes:

        building_block:
            The molecule being tested.

        binder_distance:
            The distance between binders in Angstrom.

        binder_centroid_angle:

        binder_adjacent_torsion:

        binder_angles:

        name:
            The name of the test case.

    """

    def __init__(
        self,
        building_block: stk.BuildingBlock,
        binder_distance: float,
        binder_centroid_angle: float,
        binder_adjacent_torsion: float,
        binder_angles: tuple[float, float],
        name: str,
    ) -> None:
        self.building_block = building_block
        self.binder_distance = binder_distance
        self.binder_centroid_angle = binder_centroid_angle
        self.binder_adjacent_torsion = binder_adjacent_torsion
        self.binder_angles = binder_angles
        self.name = name
