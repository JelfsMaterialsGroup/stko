"""
OpenBabel Optimizers
================

#. :class:`.OpenBabel`

Wrappers for optimizers within the `openbabel` code.

"""

import logging
import os
from openbabel import openbabel

from .optimizers import Optimizer

logger = logging.getLogger(__name__)


class OpenBabel(Optimizer):
    """
    Use OpenBabel to optimize molecules with forcefields.[1]_

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
    .. [1] http://openbabel.org/dev-api/classOpenBabel_1_1OBForceField.shtml#a2f2732698efde5c2f155bfac08fd9ded

    """

    def __init__(self, forcefield, steps=500):
        """
        Initialize `openbabel` forcefield energy calculation.

        Parameters
        ----------
        forcefield : :class:`str`
            Forcefield to use. Options include `uff`, `gaff`,
            `ghemical`, `mmff94`.

        steps : :class:`int`
            Number of steps in Conjugate Gradient optimisation.

        """

        self._forcefield = forcefield
        self._steps = steps

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
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat("mol")
        OBMol = openbabel.OBMol()
        obConversion.ReadFile(OBMol, temp_file)
        os.system('rm temp.mol')

        forcefield = openbabel.OBForceField.FindForceField(
            self._forcefield
        )
        forcefield.Setup(OBMol)
        forcefield.ConjugateGradients(self._steps)
        forcefield.GetCoordinates(OBMol)

        temp_file = 'temp.mol'
        obConversion.WriteFile(OBMol, temp_file)
        mol = mol.with_structure_from_file(temp_file)
        os.system('rm temp.mol')

        return mol
