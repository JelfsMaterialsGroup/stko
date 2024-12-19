import logging

import stk

from stko._internal.utilities.exceptions import WrapperNotInstalledError

try:
    import MDAnalysis as mda  # noqa: N813
except ModuleNotFoundError:
    mda = None

logger = logging.getLogger(__name__)


class MDAnalysis:
    """Converter for :class:`stk.Molecule` to and from MDAnalysis.

    Raises:
        :class:`WrapperNotInstalledError`: if `MDAnalysis` not installed.

    Examples:
        An stk molecule can be converted into an MDAnalysis Universe.

        .. testcode:: mdanalysis

            import stk
            import stko
            import numpy as np

            stkmol = stk.BuildingBlock('NCCNCCN').with_centroid(
                position=np.array((10, 10, 10))
            )
            universe = stko.MDAnalysis().get_universe(stkmol)


            # Can now use mdanalysis methods and analysis.
            radius_gyration = universe.atoms.radius_of_gyration()
            b_sphere = universe.atoms.bsphere()
            # These should be similar.
            universe_com = universe.atoms.center_of_mass()
            stk_centroid = stkmol.get_centroid()

        .. testcode:: mdanalysis
            :hide:

            assert np.allclose(
                universe_com, np.array([9.95245937, 10.06867664, 9.91669654])
            )
            assert np.allclose(stk_centroid, np.array([10.0, 10.0, 10.0]))

    """

    def __init__(self) -> None:
        if mda is None:
            msg = (
                "MDAnalysis is not installed; see README for " "installation."
            )
            raise WrapperNotInstalledError(msg)

    def get_universe(self, mol: stk.Molecule):  # noqa: ANN202
        """Get an MDAnalysis object.

        Parameters:
            mol:
                Molecule to convert.

        Returns:
            The MDAnalysis Universe of the molecule.

        """
        rdkit_mol = mol.to_rdkit_mol()
        return mda.Universe(rdkit_mol)
