class CaseData:
    """
    An :class:`.DummyAtom` test case.

    Attributes
    ----------
    atom : :class:`.DummyAtom`
        The atom being tested.

    id : :class:`int`
        The correct id.

    """

    def __init__(self, atom, id):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        atom : :class:`.DummyAtom`
            The atom being tested.

        id : :class:`int`
            The correct id.

        """

        self.atom = atom
        self.id = id
