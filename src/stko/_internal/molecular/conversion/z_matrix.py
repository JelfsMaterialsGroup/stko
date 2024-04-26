import logging

import stk

from stko._internal.utilities.exceptions import ConversionError
from stko._internal.utilities.utilities import (
    calculate_angle,
    calculate_dihedral,
    get_atom_distance,
)

logger = logging.getLogger(__name__)


class ZMatrix:
    """Converter for :class:`stk.Molecule` to Z-Matrix.

    Examples:
        The Z-matrix is returned as a string.

        .. code-block:: python

            import stk
            import stko

            bb1 = stk.BuildingBlock('NCCNCCN')

            print(stko.ZMatrix().get_zmatrix(bb1))

    """

    def get_zmatrix(self, molecule: stk.Molecule) -> str:
        """Get Z-matrix of a molecule.

        Parameters:
            molecule:
                Molecule to convert.

        Returns:
            The Z-matrix of the molecule.
            Distances in Angstrom, angles and torsions in degrees.

        """
        zmatrix = []
        mol_atoms = list(molecule.get_atoms())
        position_matrix = molecule.get_position_matrix()
        coords = 0
        for i, atom in enumerate(mol_atoms):
            if i == 0:
                zmatrix.append(f"{atom.__class__.__name__}")
            elif i == 1:
                distance = get_atom_distance(
                    position_matrix=position_matrix,
                    atom1_id=atom.get_id(),
                    atom2_id=mol_atoms[i - 1].get_id(),
                )

                distance = round(distance, 2)
                zmatrix.append(f"{atom.__class__.__name__} {i} {distance}")
                coords += 1
            elif i == 2:  # noqa: PLR2004
                distance = get_atom_distance(
                    position_matrix=position_matrix,
                    atom1_id=atom.get_id(),
                    atom2_id=mol_atoms[i - 1].get_id(),
                )
                angle = calculate_angle(
                    pt1=next(
                        molecule.get_atomic_positions(
                            mol_atoms[i - 3].get_id()
                        )
                    ),
                    pt2=next(
                        molecule.get_atomic_positions(
                            mol_atoms[i - 2].get_id()
                        )
                    ),
                    pt3=next(
                        molecule.get_atomic_positions(
                            mol_atoms[i - 1].get_id()
                        )
                    ),
                )

                distance = round(distance, 2)
                angle = round(angle, 2)
                zmatrix.append(
                    f"{atom.__class__.__name__} {i} {distance}"
                    f"{i-1} {angle}"
                )
                coords += 2

            else:
                distance = get_atom_distance(
                    position_matrix=position_matrix,
                    atom1_id=atom.get_id(),
                    atom2_id=mol_atoms[i - 1].get_id(),
                )
                angle = calculate_angle(
                    pt1=next(
                        molecule.get_atomic_positions(
                            mol_atoms[i - 3].get_id()
                        )
                    ),
                    pt2=next(
                        molecule.get_atomic_positions(
                            mol_atoms[i - 2].get_id()
                        )
                    ),
                    pt3=next(
                        molecule.get_atomic_positions(
                            mol_atoms[i - 1].get_id()
                        )
                    ),
                )
                torsion = calculate_dihedral(
                    pt1=next(
                        molecule.get_atomic_positions(
                            mol_atoms[i - 3].get_id()
                        )
                    ),
                    pt2=next(
                        molecule.get_atomic_positions(
                            mol_atoms[i - 2].get_id()
                        )
                    ),
                    pt3=next(
                        molecule.get_atomic_positions(
                            mol_atoms[i - 1].get_id()
                        )
                    ),
                    pt4=next(molecule.get_atomic_positions(atom.get_id())),
                )

                distance = round(distance, 2)
                angle = round(angle, 2)
                torsion = round(torsion, 2)
                zmatrix.append(
                    f"{atom.__class__.__name__} {i} {distance}"
                    f"{i-1} {angle} {i-2} {torsion}"
                )
                coords += 3

        if 3 * len(zmatrix) - 6 != coords:
            msg = (
                f"zmatrix: There are {coords}, not {3*len(zmatrix)-6} as "
                "expected. Therefore, the conversion has failed."
            )
            raise ConversionError(msg)

        return "\n".join(zmatrix)
