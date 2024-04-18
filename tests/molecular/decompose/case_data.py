import stk


class CaseData:
    """A test case.

    Attributes
    ----------
        cage:
            The cage to decompose.

        num_ligands:
            The expected number of ligands.

        metal_atom_nos:
            The metal atom atomic numbers to use.

        bb_smiles:
            The expected smiles of the extracted ligands.

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
