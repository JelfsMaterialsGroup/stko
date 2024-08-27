from copy import deepcopy

import stk

import stko
from tests.optimizers.optimizer.conftest import (
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
