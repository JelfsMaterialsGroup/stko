import logging

import numpy as np
import stk

from stko._internal.molecular.networkx.network import Network

logger = logging.getLogger(__name__)


class DecomposeMOC:
    """Decompose a metal-organic cage to obtain its organic building blocks.

    .. warning::
        This code is only present in the latest versions of stko
        that require Python 3.11!

    """

    def _get_connected_graphs(
        self,
        molecule: stk.Molecule,
        metal_atom_nos: tuple[int],
    ) -> list:
        graph = Network.init_from_molecule(molecule)
        graph = graph.with_deleted_elements(metal_atom_nos)
        return graph.get_connected_components()

    def decompose(
        self,
        molecule: stk.Molecule,
        metal_atom_nos: tuple[int],
    ) -> tuple[stk.Molecule, ...]:
        """Decompose a MOC into ligands by deleting specific metal atoms.

        Parameters:
            molecule:
                The molecule to decompose.

            metal_atom_nos:
                The atomic numbers to delete. Can be a tuple of one or
                any element on periodic table.

        Returns:
            The ligands as distinct, connected molecules.

        """
        connected_graphs = self._get_connected_graphs(
            molecule=molecule,
            metal_atom_nos=metal_atom_nos,
        )
        ligands = []
        for cg in connected_graphs:
            # Get atoms from nodes.
            atoms = list(cg)
            atom_ids = tuple(i.get_id() for i in atoms)

            # Sort both by atom id.
            atom_ids, atoms = zip(  # type: ignore[assignment]
                *sorted(zip(atom_ids, atoms, strict=True)), strict=True
            )

            atom_ids_map = {atom_ids[i]: i for i in range(len(atom_ids))}
            new_mol = stk.BuildingBlock.init(
                atoms=(
                    stk.Atom(
                        id=atom_ids_map[i.get_id()],
                        atomic_number=i.get_atomic_number(),
                        charge=i.get_charge(),
                    )
                    for i in atoms
                ),
                bonds=(
                    i.with_ids(id_map=atom_ids_map)
                    for i in molecule.get_bonds()
                    if i.get_atom1().get_id() in atom_ids
                    and i.get_atom2().get_id() in atom_ids
                ),
                position_matrix=np.array(
                    tuple(
                        i
                        for i in molecule.get_atomic_positions(
                            atom_ids=atom_ids
                        )
                    )
                ),
            )

            ligands.append(new_mol)

        return tuple(ligands)
