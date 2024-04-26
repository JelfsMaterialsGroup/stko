from dataclasses import dataclass

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    building_block: stk.BuildingBlock
    binder_distance: float
    binder_centroid_angle: float
    binder_adjacent_torsion: float
    binder_angles: tuple[float, float]
    binder_binder_angle: float
    name: str

    @property
    def bite_angles(self) -> tuple[float, float]:
        return self.binder_angles[0] - 90, self.binder_angles[1] - 90
