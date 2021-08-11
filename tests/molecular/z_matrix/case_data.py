
class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`stk.Molecule`
        The molecule to test.

    zmatrix : :class:`str`
        The correct Z-matrix string.

    """

    def __init__(self, molecule, zmatrix):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`stk.Molecule`
            The molecule to test.

        zmatrix : :class:`str`
            The correct Z-matrix string.

        """

        self.molecule = molecule
        self.zmatrix = zmatrix
