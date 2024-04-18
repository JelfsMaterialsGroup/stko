class CaseData:
    """An :class:`.Du` test case.

    Attributes:
    ----------
    atom : :class:`.Du`
        The atom being tested.

    id : :class:`int`
        The correct id.

    """

    def __init__(self, atom, id):
        """Parameters
        ----------
        atom : :class:`.Du`
            The atom being tested.

        id : :class:`int`
            The correct id.

        """
        self.atom = atom
        self.id = id
