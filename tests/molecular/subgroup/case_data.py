import stk


class CaseData:
    """A test case.

    Attributes
    ----------
        building_block:
            The molecule being tested.

        binder_distance:
            The distance between binders in Angstrom.

        name:
            The name of the test case.

    """

    def __init__(
        self,
        molecule: stk.Molecule,
        sub_group_data: dict[str, float],
        name: str,
    ) -> None:
        self.molecule = molecule
        self.sub_group_data = sub_group_data
        self.name = name
