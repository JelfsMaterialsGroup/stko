import stk


class CaseData:
    """
    A test case.

    Attributes:

        constructed_molecule:
            The molecule being tested.

        name:
            The name of the test case.

    """

    def __init__(
        self,
        constructed_molecule: stk.ConstructedMolecule,
        name: str,
    ) -> None:
        self.constructed_molecule = constructed_molecule
        self.name = name
