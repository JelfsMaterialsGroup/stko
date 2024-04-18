import logging

import stk
from stko._internal.utilities.exceptions import WrapperNotInstalledError

try:
    import MDAnalysis as mda
except ModuleNotFoundError:
    mda = None

logger = logging.getLogger(__name__)


class MDAnalysis:
    """Converter for :class:`stk.Molecule` to and from MDAnalysis.

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

    def __init__(self) -> None:
        """Raises
        :class:`WrapperNotInstalledError` if `MDAnalysis` not installed.


        """
        if mda is None:
            raise WrapperNotInstalledError(
                "MDAnalysis is not installed; see README for " "installation."
            )

    def get_universe(self, mol: stk.Molecule):  # type: ignore[no-untyped-def]
        """Get an MDAnalysis object.

        Parameters
        ----------
            mol:
                Molecule to convert.

        Returns
        -------
            :class:`MDAnalysis.Universe`:
                The MDAnalysis Universe of the molecule.

        """
        rdkit_mol = mol.to_rdkit_mol()
        return mda.Universe(rdkit_mol)
