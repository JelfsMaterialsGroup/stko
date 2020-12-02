# This is an example for incorporating the benzene fixture
import numpy as np
from .utilities import compare_molecules

import stko


def test_fake_opt(benzene_build):

    # Perform optimisation.
    optimizer = stko.ETKDG()
    opt_benzene = optimizer.optimize(benzene_build)

    compare_molecules(
        initial_molecule=benzene_build,
        optimized_molecule=opt_benzene,
    )
