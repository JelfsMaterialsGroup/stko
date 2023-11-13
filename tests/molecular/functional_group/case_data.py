import stko


class CaseData:
    """
    A test case.

    Attributes:

        smiles:
            The SMILES string of molecule being tested.

        num_fgs:
            The expected number of functional groups.

        factory:
            The factory to use.

        name:
            The name of the test case.

    """

    def __init__(
        self,
        smiles: str,
        num_fgs: int,
        factory: stko.functional_groups.ThreeSiteFactory,
        name: str,
    ) -> None:
        self.smiles = smiles
        self.num_fgs = num_fgs
        self.factory = factory
        self.name = name
