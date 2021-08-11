"""
Molecule Modifier
=================

#. :class:`.MoleculeModifier`

Base class for modifying a molecule.

"""

import logging

logger = logging.getLogger(__name__)


class MoleculeModifier:
    """
    Modify a molecule into one or many new molecules.

    """

    def modify(self, molecule):
        """
        Modify a molecule.

        Parameters
        ----------
        molecule : :class:`stk.Molecule`
            Molecule to modify.

        Returns
        -------
        molecules : :class:`iterable` of :class:`stk.BuildingBlock`
            The resulting list of molecules.

        """

        raise NotImplementedError()
