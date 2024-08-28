import pytest
import stk
from pytest_lazyfixture import lazy_fixture

import stko

from .case_data import CaseData


@pytest.fixture
def case_data_1(
    atomic_number: int,
    id: int,  # noqa: A002
    charge: int,
    position: tuple[float, float, float],
) -> CaseData:
    """A :class:`.CaseData` instance."""
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


@pytest.fixture(params=(lazy_fixture("case_data_1"),))
def case_data(request: pytest.FixtureRequest) -> CaseData:
    """A :class:`.CaseData` instance."""
    return request.param


@pytest.fixture(
    params=[0, 3],
)
def id(request: pytest.FixtureRequest) -> int:  # noqa: A001
    """An atom id."""
    return request.param


@pytest.fixture(
    params=[0],
)
def charge(request: pytest.FixtureRequest) -> int:
    """An atomic charge."""
    return request.param


@pytest.fixture(params=tuple(range(1, 118)))
def atomic_number(request: pytest.FixtureRequest) -> int:
    """An atomic number."""
    return request.param


@pytest.fixture(
    params=[(0, 0, 0), (100, 0, -1)],
)
def position(request: pytest.FixtureRequest) -> tuple[float, float, float]:
    """A position."""
    return request.param


@pytest.fixture(
    params=[
        cls
        for cls in stk.__dict__.values()
        if isinstance(cls, type)
        and issubclass(cls, stk._internal.elements.AtomImpl)  # noqa: SLF001
        and cls is not stk._internal.elements.AtomImpl  # noqa: SLF001
    ],
)
def cls(request: pytest.FixtureRequest):  # noqa: ANN201
    return request.param


@pytest.fixture
def positioned_atom(
    cls,  # noqa: ANN001
    position: tuple[float, float, float],
) -> stko.PositionedAtom:
    """An :class:`.PositionedAtom` instance."""
    return stko.PositionedAtom(
        atom=cls(3, -5),
        position=position,
    ).clone()
