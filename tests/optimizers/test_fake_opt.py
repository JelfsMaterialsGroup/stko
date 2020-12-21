from .utilities import compare_benzenes

import stko


def test_fake_opt(benzene_build):

    # Perform optimisation.
    optimizer = stko.ETKDG()
    opt_benzene = optimizer.optimize(benzene_build)

    compare_benzenes(
        initial_molecule=benzene_build,
        optimized_molecule=opt_benzene,
    )
