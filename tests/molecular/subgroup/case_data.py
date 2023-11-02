import stk


class CaseData:
    """
    A test case.

    Attributes:

        building_block:
            The molecule being tested.

        binder_distance:
            The distance between binders in Angstrom.

        name:
            The name of the test case.

    """

    def __init__(
        self,
        building_block: stk.BuildingBlock,
        binder_distance: float,
        name: str,
    ) -> None:
        self.building_block = building_block
        self.binder_distance = binder_distance
        self.name = name
