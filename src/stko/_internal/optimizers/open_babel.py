import logging
from pathlib import Path

import numpy as np

try:
    from openbabel import openbabel
except ImportError:
    openbabel = None

from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.types import MoleculeT
from stko._internal.utilities.exceptions import (
    ForceFieldSetupError,
    WrapperNotInstalledError,
)

logger = logging.getLogger(__name__)


class OpenBabel(Optimizer):
    """Use OpenBabel to optimize molecules with forcefields.

    Parameters:
        forcefield:
            Forcefield to use. Options include `uff`, `gaff`,
            `ghemical`, `mmff94`.

        repeat_steps:
            Number of optimisation steps. Each optimisation step
            contains `sd_steps` steepest descent and then `cg_steps`
            conjugate gradient runs.

        sd_steps:
            Number of steepest descent steps per optimisations.

        cg_steps:
                Number of conjugate gradient steps per optimisations.

    Raises:
        :class:`WrapperNotInstalledError`:
            if `openbabel` not installed.

    .. warning::
        this optimizer seems to be machine dependant, producing
        different energies after optimisation on Ubunut 18 vs. Ubuntu 20.

    See Also:
        * OpenBabel: https://github.com/openbabel/openbabel

    Examples:
        .. code-block:: python

            import stk
            import stko

            mol = stk.BuildingBlock('NCCNCCN')
            openbabel = stko.OpenBabel('uff')
            mol = openbabel.optimize(mol)

    """

    def __init__(
        self,
        forcefield: str,
        repeat_steps: int = 10,
        sd_steps: int = 50,
        cg_steps: int = 50,
    ) -> None:
        if openbabel is None:
            msg = "openbabel is not installed; see README for " "installation."
            raise WrapperNotInstalledError(msg)

        self._forcefield = forcefield
        self._repeat_steps = repeat_steps
        self._sd_steps = sd_steps
        self._cg_steps = cg_steps

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        temp_file = Path("temp.mol")
        mol.write(temp_file)
        try:
            ob_conversion = openbabel.OBConversion()
            ob_conversion.SetInFormat("mol")
            ob_mol = openbabel.OBMol()
            ob_conversion.ReadFile(ob_mol, temp_file)
            ob_mol.PerceiveBondOrders()
        finally:
            temp_file.unlink()

        forcefield = openbabel.OBForceField.FindForceField(self._forcefield)
        outcome = forcefield.Setup(ob_mol)
        if not outcome:
            msg = f"Openbabel: {self._forcefield} could not be setup for {mol}"
            raise ForceFieldSetupError(msg)
        for _ in range(self._repeat_steps):
            forcefield.SteepestDescent(self._sd_steps)
            forcefield.GetCoordinates(ob_mol)
            forcefield.ConjugateGradients(self._cg_steps)
            forcefield.GetCoordinates(ob_mol)

        position_matrix = [
            np.array([atom.GetX(), atom.GetY(), atom.GetZ()])
            for atom in openbabel.OBMolAtomIter(ob_mol)
        ]

        return mol.with_position_matrix(np.asarray(position_matrix))
