import stk


class CaseData:
    def __init__(
        self,
        molecule: stk.Molecule,
        name: str,
    ) -> None:
        self.molecule = molecule
        self.name = name
