class CaseData:
    """
    An :class:`.PositionedAtom` test case.

    Attributes
    ----------
    atom : :class:`.PositionedAtom`
        The atom being tested.

    id : :class:`int`
        The correct id.

    charge : :class:`int`
        The correct charge.

    atomic_number : :class:`int`
        The correct atomic number.

    position : :class:`tuple` of :class:`float`
        The position (`x`, `y`, `z`) of the atom in cartesian
        coordinates.

    """

    def __init__(self, atom, id, charge, atomic_number, position):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        atom : :class:`.PositionedAtom`
            The atom being tested.

        id : :class:`int`
            The correct id.

        charge : :class:`int`
            The correct charge.

        atomic_number : :class:`int`
            The correct atomic number.

        position : :class:`tuple` of :class:`float`
            The position (`x`, `y`, `z`) of the atom in cartesian
            coordinates.

        """

        self.atom = atom
        self.id = id
        self.charge = charge
        self.atomic_number = atomic_number
        self.position = position
