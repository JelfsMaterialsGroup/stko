"""
OpenBabel Optimizers
================

#. :class:`.OpenBabel`

Wrappers for optimizers within the `openbabel` code.

"""

import logging
import os
import numpy as np
try:
    from openbabel import openbabel
except ImportError:
    openbabel = None

from .optimizers import Optimizer
from ..utilities import WrapperNotInstalledException

logger = logging.getLogger(__name__)


class OpenBabelError(Exception):
    ...


class ForceFieldSetupError(OpenBabelError):
    ...


class OpenBabel(Optimizer):
    """
    Use OpenBabel to optimize molecules with forcefields.[1]_

    Warning: this optimizer seems to be machine dependant, producing
    different energies after optimisation on Ubunut 18 vs. Ubuntu 20.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        mol = stk.BuildingBlock('NCCNCCN')
        openbabel = stko.OpenBabel('uff')
        mol = openbabel.optimize(mol)

    References
    ----------
    .. [1] http://openbabel.org/dev-api/classOpenBabel_1_
    1OBForceField.shtml#a2f2732698efde5c2f155bfac08fd9ded

    """

    def __init__(
        self,
        forcefield,
        repeat_steps=10,
        sd_steps=50,
        cg_steps=50,
    ):
        """
        Initialize `openbabel` forcefield energy calculation.

        Parameters
        ----------
        forcefield : :class:`str`
            Forcefield to use. Options include `uff`, `gaff`,
            `ghemical`, `mmff94`.

        repeat_steps : :class:`int`
            Number of optimisation steps. Each optimisation step
            contains `sd_steps` steepest descent and then `cg_steps`
            conjugate gradient runs.

        sd_steps : :class:`int`
            Number of steepest descent steps per optimisations.

        cg_steps : :class:`int`
            Number of conjugate gradient steps per optimisations.

        """

        if openbabel is None:
            raise WrapperNotInstalledException(
                'openbabel is not installed; see README for '
                'installation.'
            )

        self._forcefield = forcefield
        self._repeat_steps = repeat_steps
        self._sd_steps = sd_steps
        self._cg_steps = cg_steps

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        mol : :class:`.Molecule`
            The optimized molecule.

        """

        temp_file = 'temp.mol'
        mol.write(temp_file)
        try:
            obConversion = openbabel.OBConversion()
            obConversion.SetInFormat("mol")
            OBMol = openbabel.OBMol()
            obConversion.ReadFile(OBMol, temp_file)
            OBMol.PerceiveBondOrders()
        finally:
            os.system('rm temp.mol')

        forcefield = openbabel.OBForceField.FindForceField(
            self._forcefield
        )
        outcome = forcefield.Setup(OBMol)
        if not outcome:
            raise ForceFieldSetupError(
                f"{self._forcefield} could not be setup for {mol}"
            )
        for step in range(self._repeat_steps):
            forcefield.SteepestDescent(self._sd_steps)
            forcefield.GetCoordinates(OBMol)
            forcefield.ConjugateGradients(self._cg_steps)
            forcefield.GetCoordinates(OBMol)

        position_matrix = []
        for atom in openbabel.OBMolAtomIter(OBMol):
            # get coordinates
            position_matrix.append(
                np.array([atom.GetX(), atom.GetY(), atom.GetZ()])
            )

        mol = mol.with_position_matrix(np.asarray(position_matrix))
        return mol
