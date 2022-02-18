"""
MD Analysis
===========

#. :class:`.MDAnalysis`

Class for converting a molecule to and back from an MDAnalysis object.

"""

import logging

import MDAnalysis as mda

logger = logging.getLogger(__name__)


class MDAnalysis:
    """
    Converter for :class:`stk.Molecule` to and from MDAnalysis.

    Examples
    --------

    An stk molecule can be converted into an MDAnalysis Universe.

    .. code-block:: python

        import stk
        import stko

        stkmol = stk.BuildingBlock('NCCNCCN').with_centroid(
            position=np.array((10, 10, 10))
        )
        universe = stko.MDAnalysis().get_universe(stkmol)

        print('R_g:', universe.atoms.radius_of_gyration())
        print('B_sphere:', universe.atoms.bsphere())
        print('Universe COM:', universe.atoms.center_of_mass())
        print('stk centroid:', stkmol.get_centroid())

    """

    def get_universe(self, mol):
        """
        Get an MDAnalysis object.

        Parameters
        ----------
        mol : :class:`stk.Molecule`
            Molecule to convert.

        Returns
        -------
        :class:`MDAnalysis.Universe`
            The MDAnalysis Universe of the molecule.

        """

        rdkit_mol = mol.to_rdkit_mol()
        return mda.Universe(rdkit_mol)
