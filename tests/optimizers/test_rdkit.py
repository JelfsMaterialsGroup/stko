import rdkit.Chem as Chem
from stko.optimizers import MMFF, UFF, ETKDG
import stk

from ..fixtures.benzene import benzene_build
import numpy as np
from .utilities import compare_molecules

import stko

def test_MMFF_opt(benzene_build):

    # Perform optimisation.
    optimizer = stko.MMFF()
    opt_benzene = optimizer.optimize(benzene_build)

    compare_molecules(
        initial_molecule=benzene_build,
        optimized_molecule=opt_benzene,
    )

def test_UFF_opt(benzene_build):

    # Perform optimisation.
    optimizer = stko.UFF()
    opt_benzene = optimizer.optimize(benzene_build)

    compare_molecules(
        initial_molecule=benzene_build,
        optimized_molecule=opt_benzene,
    )

def test_ETKDG_opt(benzene_build):

    # Perform optimisation.
    optimizer = stko.ETKDG()
    opt_benzene = optimizer.optimize(benzene_build)

    compare_molecules(
        initial_molecule=benzene_build,
        optimized_molecule=opt_benzene,
    )
