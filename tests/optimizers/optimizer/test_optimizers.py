import pathlib
from copy import deepcopy

import pytest
import stk

import stko
from tests.optimizers.optimizer.conftest import (
    DummyFileIOOptimizer,
    FailingOptimizer,
    PassingOptimizer,
)
from tests.optimizers.utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)


def test_optimizer_sequence(
    passing_optimizer: PassingOptimizer,
    unoptimized_mol: stk.BuildingBlock,
) -> None:
    opts = [deepcopy(passing_optimizer) for _ in range(10)]
    opt_seq = stko.OptimizerSequence(*opts)
    opt_res = opt_seq.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)


def test_trycatchoptimizer(
    passing_optimizer: PassingOptimizer,
    failing_optimizer: FailingOptimizer,
    unoptimized_mol: stk.BuildingBlock,
) -> None:
    opt = stko.TryCatchOptimizer(
        try_optimizer=failing_optimizer,
        catch_optimizer=passing_optimizer,
    )
    opt_res = opt.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)


def test_fileiooptimizer_none_output_dir_creates_directory(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    opt = DummyFileIOOptimizer(output_dir=None)
    out = opt._setup_output_dir()  # noqa: SLF001
    assert out.is_dir()


def test_fileiooptimizer_existing_output_dir_raises_fileexistserror(
    tmp_path: pathlib.Path,
) -> None:
    output_dir = tmp_path / "test_output_dir"
    output_dir.mkdir(parents=True)

    opt = DummyFileIOOptimizer(
        output_dir=output_dir,
        delete_path=False,
    )
    with pytest.raises(FileExistsError):
        opt._setup_output_dir()  # noqa: SLF001


def test_fileiooptimizer_cwd_output_dir_raises_valueerror() -> None:
    output_dir = pathlib.Path.cwd()

    with pytest.raises(
        ValueError,
        match="Output directory cannot be the current working directory",
    ):
        DummyFileIOOptimizer(
            output_dir=output_dir,
        )
