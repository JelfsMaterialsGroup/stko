import logging
from pathlib import Path
from typing import Self

import numpy as np
from stk import PeriodicInfo

from stko._internal.molecular.periodic.utilities import get_from_parameters

logger = logging.getLogger(__name__)


class UnitCell(PeriodicInfo):
    """Unit cell information for periodic systems.

    We are aware that this naming choice may not be appropriate (because
    not all inputs will be unit cells, strictly). However, for backwards
    compatability, we have not changed this naming.

    """

    @classmethod
    def _update_periodic_info(
        cls,
        vector_1: np.ndarray,
        vector_2: np.ndarray,
        vector_3: np.ndarray,
    ) -> Self:
        """Return clone of :class:`.UnitCell` with new parameters."""
        clone = cls.__new__(cls)
        UnitCell.__init__(
            self=clone,
            vector_1=vector_1,
            vector_2=vector_2,
            vector_3=vector_3,
        )

        return clone

    def with_cell_from_vectors(
        self,
        vector_1: np.ndarray,
        vector_2: np.ndarray,
        vector_3: np.ndarray,
    ) -> Self:
        """Update cell.

        Parameters:
            vector_1:
                First cell lattice vector of shape (3, ) in
                Angstrom.

            vector_2:
                Second cell lattice vector of shape (3, ) in
                Angstrom.

            vector_3:
                Third cell lattice vector of shape (3, ) in
                Angstrom.

        Returns:
            Clone with updated cell parameters.

        """
        return self._update_periodic_info(
            vector_1=vector_1,
            vector_2=vector_2,
            vector_3=vector_3,
        )

    def with_cell_from_turbomole(self, filename: Path | str) -> Self:  # noqa:  PLR0912, C901
        """Update cell from structure in Turbomole coord file.

        Returns:
            Clone with updated cell parameters.

        """
        filename = Path(filename)
        bohr_to_ang = 0.5291772105638411

        with filename.open() as f:
            content = f.readlines()

        periodicity: bool | int = False
        lattice_vectors = None
        lattice_units = None
        cell_parameters = None
        cell_units = None
        for line_number, line in enumerate(content):
            if "$periodic" in line:
                periodicity = int(line.rstrip().split()[1])
            if "$cell" in line:
                if "angs" in line:
                    cell_units = "angstrom"
                elif "bohr" in line:
                    cell_units = "bohr"
                else:
                    msg = "cell not in Angstroms."
                    raise ValueError(msg)
                cell_parameters = [
                    float(j) for j in content[line_number + 1].rstrip().split()
                ]
            if "$lattice" in line:
                if "angs" in line:
                    lattice_units = "angstrom"
                elif "bohr" in line:
                    lattice_units = "bohr"
                else:
                    msg = "lattice not in Angstroms."
                    raise ValueError(msg)
                lattice_vectors = (
                    np.array(
                        [
                            float(j)
                            for j in content[line_number + 1].rstrip().split()
                        ]
                    ),
                    np.array(
                        [
                            float(j)
                            for j in content[line_number + 2].rstrip().split()
                        ]
                    ),
                    np.array(
                        [
                            float(j)
                            for j in content[line_number + 3].rstrip().split()
                        ]
                    ),
                )

        # Check that cell is only defined once.
        chk2 = lattice_vectors is not None and cell_parameters is not None
        if periodicity and chk2:
            msg = "The cell is defined twice in the file."
            raise RuntimeError(msg)

        if lattice_vectors is not None:
            vector_1 = (
                lattice_vectors[0] * bohr_to_ang
                if lattice_units == "bohr"
                else lattice_vectors[0]
            )
            vector_2 = (
                lattice_vectors[1] * bohr_to_ang
                if lattice_units == "bohr"
                else lattice_vectors[0]
            )
            vector_3 = (
                lattice_vectors[2] * bohr_to_ang
                if lattice_units == "bohr"
                else lattice_vectors[0]
            )
        elif cell_parameters is not None:
            vector_1, vector_2, vector_3 = get_from_parameters(
                a=(
                    cell_parameters[0] * bohr_to_ang
                    if cell_units == "bohr"
                    else cell_parameters[0]
                ),
                b=(
                    cell_parameters[1] * bohr_to_ang
                    if cell_units == "bohr"
                    else cell_parameters[1]
                ),
                c=(
                    cell_parameters[2] * bohr_to_ang
                    if cell_units == "bohr"
                    else cell_parameters[2]
                ),
                alpha=cell_parameters[3],
                beta=cell_parameters[4],
                gamma=cell_parameters[5],
            )
        else:
            msg = "The cell is not defined in the file."
            raise RuntimeError(msg)
        # Update the cell.
        return self._update_periodic_info(vector_1, vector_2, vector_3)

    def with_cell_from_cif(self, filename: Path | str) -> Self:
        """Update cell from structure in CIF.

        Returns:
            Clone with updated cell parameters.

        """
        filename = Path(filename)
        cell_info: dict[str, float] = {}

        targets = {
            "_cell_length_a": "a",
            "_cell_length_b": "b",
            "_cell_length_c": "c",
            "_cell_angle_alpha": "alpha",
            "_cell_angle_beta": "beta",
            "_cell_angle_gamma": "gamma",
        }

        with filename.open() as f:
            lines = f.readlines()

        for targ in targets:
            for line in lines:
                # Avoid running through the rest.
                if targets[targ] in cell_info:
                    break
                splits = line.rstrip().split(" ")
                if splits[0] == targ:
                    cell_info[targets[targ]] = float(splits[-1])

        vector_1, vector_2, vector_3 = get_from_parameters(
            a=cell_info["a"],
            b=cell_info["b"],
            c=cell_info["c"],
            alpha=cell_info["alpha"],
            beta=cell_info["beta"],
            gamma=cell_info["gamma"],
        )
        # Update the cell.
        return self._update_periodic_info(vector_1, vector_2, vector_3)
