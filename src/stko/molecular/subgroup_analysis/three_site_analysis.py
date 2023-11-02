import logging

import stk

from stko.functional_groups import ThreeSiteFG
from stko.utilities.exceptions import NotDitopicThreeSiteError
from stko.utilities.utilities import (
    get_atom_distance,
)

logger = logging.getLogger(__name__)


class DitopicThreeSiteAnalyser:
    """
    Analyses geometry of functional groups in ditopic molecules.

    """

    def _check_functional_groups(self, molecule: stk.BuildingBlock) -> None:
        """
        Check if the molecule has two ditopic functional groups.

        Raises:

            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """

        fg_counts = 0
        for fg in molecule.get_functional_groups():
            if isinstance(fg, ThreeSiteFG):
                fg_counts += 1

        if fg_counts != 2:
            raise NotDitopicThreeSiteError(
                f"{molecule} does not have 2 ThreeSiteFG functional groups."
            )

    def get_binder_distance(self, molecule: stk.BuildingBlock) -> float:
        """
        Get the distance between binder atoms in Angstrom.

        Raises:

            NotDitopicThreeSiteError: if does not have two ThreeSiteFG.

        """

        self._check_functional_groups(molecule)

        binder_ids = [
            fg.get_binder().get_id()  # type: ignore[attr-defined]
            for fg in molecule.get_functional_groups()
        ]
        return get_atom_distance(
            position_matrix=molecule.get_position_matrix(),
            atom1_id=binder_ids[0],
            atom2_id=binder_ids[1],
        )
