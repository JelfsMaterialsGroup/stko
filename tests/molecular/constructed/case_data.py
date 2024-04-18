import numpy as np
import stk


class CaseData:
    """A test case.

    Attributes
    ----------
        constructed_molecule:
            The molecule being tested.

        name:
            The name of the test case.

    """

    def __init__(
        self,
        constructed_molecule: stk.ConstructedMolecule,
        centroids: dict[int, np.ndarray],
        atom_ids: dict[int, list[int]],
        name: str,
    ) -> None:
        self.constructed_molecule = constructed_molecule
        self.centroids = centroids
        self.atom_ids = atom_ids
        self.name = name
