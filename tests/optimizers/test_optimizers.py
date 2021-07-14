import stko
from copy import deepcopy
from .utilities import (
    inequivalent_position_matrices,
    is_equivalent_molecule,
)


def test_optimizer_sequence(
    passing_optimizer, unoptimized_mol
):

    opts = [deepcopy(passing_optimizer) for i in range(10)]
    opt_seq = stko.OptimizerSequence(*opts)
    opt_res = opt_seq.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)


def test_trycatchoptimizer(
    passing_optimizer, failing_optimizer, unoptimized_mol,
):

    opt = stko.TryCatchOptimizer(
        try_optimizer=failing_optimizer,
        catch_optimizer=passing_optimizer,
    )
    opt_res = opt.optimize(unoptimized_mol)
    is_equivalent_molecule(opt_res, unoptimized_mol)
    inequivalent_position_matrices(opt_res, unoptimized_mol)
