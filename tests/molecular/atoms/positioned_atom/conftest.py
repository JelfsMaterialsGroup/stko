import pytest
from pytest_lazyfixture import lazy_fixture
import stk
import stko

from .case_data import CaseData


@pytest.fixture
def case_data_1(atomic_number, id, charge, position):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        atom=stko.PositionedAtom(
            atom=stk.Atom(
                id=id,
                atomic_number=atomic_number,
                charge=charge,
            ),
            position=position,
        ),
        id=id,
        charge=charge,
        atomic_number=atomic_number,
        position=position,
    )


@pytest.fixture(
    params=(
        lazy_fixture('case_data_1'),
    )
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param


@pytest.fixture(
    params=[0, 3],
)
def id(request):
    """
    An atom id.

    """

    return request.param


@pytest.fixture(
    params=[0],
)
def charge(request):
    """
    An atomic charge.

    """

    return request.param


@pytest.fixture(params=tuple(range(1, 118)))
def atomic_number(request):
    """
    An atomic number.

    """

    return request.param


@pytest.fixture(
    params=[(0, 0, 0), (100, 0, -1)],
)
def position(request):
    """
    A position.

    """

    return request.param


@pytest.fixture(
    params=[
        cls for cls in stk.__dict__.values()
        if isinstance(cls, type)
        and issubclass(cls, stk.AtomImpl)
        and cls is not stk.AtomImpl
    ],
)
def cls(request):
    """
    Return an :class:`.Atom` instance.

    Parameters
    ----------
    id : :class:`int`
        The id of the returned atom.

    charge : :class:`float`
        The charge of the returned atom.

    Returns
    -------
    :class:`.Atom`
        An atom.

    """

    return request.param


@pytest.fixture
def positioned_atom(cls, position):
    """
    An :class:`.PositionedAtom` instance.

    """

    return stko.PositionedAtom(
        atom=cls(3, -5),
        position=position,
    ).clone()
