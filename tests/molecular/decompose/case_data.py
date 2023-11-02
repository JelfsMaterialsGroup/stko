import stk


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
        cage: stk.Molecule,
        num_ligands: int,
        metal_atom_nos: tuple[int],
        bb_smiles: tuple[str],
        name: str,
    ) -> None:
        self.cage = cage
        self.num_ligands = num_ligands
        self.metal_atom_nos = metal_atom_nos
        self.bb_smiles = bb_smiles
        self.name = name
