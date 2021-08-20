import pytest
from pytest_lazyfixture import lazy_fixture
import stko

from .case_data import CaseData


@pytest.fixture
def case_data_1(id):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        atom=stko.Du(id),
        id=id,
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


@pytest.fixture
def dummy_atom(id):
    """
    An :class:`.Du` instance.

    """

    return stko.Du(id).clone()
