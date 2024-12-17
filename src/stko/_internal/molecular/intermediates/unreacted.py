import itertools as it
import logging
from collections import abc
from dataclasses import dataclass
from functools import partial

import numpy as np
import stk

from stko._internal.molecular.networkx.network import Network

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class NamedIntermediate:
    """Named intermediate container."""

    intermediate_name: str
    molecule: stk.BuildingBlock
    present_bbs: dict[stk.BuildingBlock, int]
    num_atoms: int
    count_bbs: int

    def __str__(self) -> str:
        """String representation."""
        return repr(self)

    def __repr__(self) -> str:
        """String representation."""
        return f"{self.__class__.__name__}({self.intermediate_name})"


class UnreactedTopologyGraph:
    """A class containing topology graphs and performing subset of reactions.

    Use this to get partially reacted topology graphs.

    .. testcode:: unreacted-topology-graph

        import stk
        import stko
        import pathlib

        smarts = '[#6]~[#7]'

        mol = stk.BuildingBlock(
            smiles='C1=CC=C(C=C1)CN',
            # Find the functional groups uses to split the molecule.
            functional_groups=(stk.SmartsFunctionalGroupFactory(
                smarts=smarts,
                bonders=(),
                deleters=(),
            ), )
        )

        # For each functional group, we add disconnecting atom ids to two
        # objects for passing to extractor.
        broken_bonds_by_id = []
        disconnectors = []
        for fg in mol.get_functional_groups():
            atom_ids = list(fg.get_atom_ids())
            bond_c = atom_ids[0]
            bond_n = atom_ids[1]
            broken_bonds_by_id.append(sorted((bond_c, bond_n)))
            disconnectors.extend((bond_c, bond_n))

        new_topology_graph = stko.TopologyExtractor()
        tg_info = new_topology_graph.extract_topology(
            molecule=mol,
            broken_bonds_by_id=broken_bonds_by_id,
            disconnectors=set(disconnectors),
        )

        # Use these variables to write new topology graph like seen in
        # stk source code.
        vertex_positions = tg_info.get_vertex_positions()
        connectivities = tg_info.get_connectivities()
        edge_pairs = tg_info.get_edge_pairs()

        # Write molecules to file to visualise new topology graph.
        output_directory = pathlib.Path("output_directory")
        output_directory.mkdir(exist_ok=True)

        # Original molecule.
        mol.write(output_directory / "tg_cage.mol")
        # Underlying topology graph.
        tg_info.write(output_directory / "tg_info.pdb")
    .. testcode:: analysing-cage
        :hide:

        assert neigh_counts[0][0][3] == 17
    .. testcleanup:: unreacted-topology-graph

        import shutil

        shutil.rmtree('output_directory')

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles="C1CCC(C(C1)N)N",
            functional_groups=[stk.PrimaryAminoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles="C1=C(C=C(C=C1C=O)C=O)C=O",
            functional_groups=[stk.AldehydeFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cage.get_atoms(),
                    cage.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in cage.get_bonds()
            ),
        )

    """

    def __init__(self, topology_graph: stk.TopologyGraph) -> None:
        """Initialize UnreactedTopologyGraph.

        Parameters:
            topology_graph:
                The `stk` topology graph to contain.

        """
        self._topology_graph = topology_graph
        self._state = self._topology_graph._get_construction_state()  # noqa: SLF001
        self._state = self._topology_graph._place_building_blocks(self._state)  # noqa: SLF001
        get_reaction = partial(
            self._topology_graph._reaction_factory.get_reaction,  # noqa: SLF001
            self._state,
        )
        self._reactions = tuple(
            map(get_reaction, self._topology_graph._edge_groups)  # noqa: SLF001
        )
        self._results = tuple(
            reaction.get_result() for reaction in self._reactions
        )

    def get_available_reactions(self) -> abc.Sequence[stk.Reaction]:
        """Get all the reaction classes possible."""
        return self._reactions

    def get_reaction_results(self) -> abc.Sequence[stk.ReactionResult]:
        """Get all the reaction results."""
        return self._results

    def yield_constructed_molecules(
        self,
        n: int | None = None,
    ) -> abc.Sequence[stk.ConstructedMolecule]:
        """Yield constructed molecules for possible reaction combos.

        If `n` is None, this produces all reactions, which could be
        a combinatorial nightmare.

        """
        for i in range(1, len(self._results)):
            if n is not None and i != n:
                continue
            for inter_reactions in it.combinations(self._reactions, i):
                results = tuple(
                    reaction.get_result() for reaction in inter_reactions
                )

                yield stk.ConstructedMolecule.init_from_construction_result(
                    self._topology_graph._get_construction_result(  # noqa: SLF001
                        self._state.with_reaction_results(
                            inter_reactions, results
                        )
                    )
                )

    def separate_molecule(
        self, molecule: stk.Molecule
    ) -> abc.Sequence[stk.Molecule]:
        """Given a molecule, it returns distinct disconnected molecules."""
        network = Network.init_from_molecule(molecule)
        connected = network.get_connected_components()
        molecules = []
        for cg in connected:
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
            molecules.append(new_mol)
        return molecules

    def get_reacted_smiles(self, n: int | None = None) -> set[str]:
        """Yield constructed molecules with n reactions performed."""
        yielded_smiles = set()
        for const_mol in self.yield_constructed_molecules(n=n):
            distinct_molecules = self.separate_molecule(const_mol)
            for dmol in distinct_molecules:
                smiles = stk.Smiles().get_key(dmol)

                if "." in smiles:
                    msg = "Found `.` in smiles."
                    raise RuntimeError(msg)

                if smiles in yielded_smiles:
                    continue
                yielded_smiles.add(smiles)

        return yielded_smiles

    def get_present_building_blocks(
        self,
        const_mol: stk.ConstructedMolecule,
        subset: stk.Molecule,
    ) -> dict[stk.BuildingBlock, int]:
        """Get the building blocks present in a constructed molecule."""
        bbs = {}
        for atom_info in const_mol.get_atom_infos():
            if atom_info.get_atom().get_id() in tuple(
                i.get_id() for i in subset.get_atoms()
            ):
                if atom_info.get_building_block() not in bbs:
                    bbs[atom_info.get_building_block()] = []
                if (
                    atom_info.get_building_block_id()
                    not in bbs[atom_info.get_building_block()]
                ):
                    bbs[atom_info.get_building_block()].append(
                        atom_info.get_building_block_id()
                    )
        return bbs

    def get_named_intermediates(self, n: int | None = None) -> set[str]:
        """Yield constructed molecules with up to n reactions performed."""
        yielded_smiles = set()
        intermediates = []
        for const_mol in self.yield_constructed_molecules(n=n):
            distinct_molecules = self.separate_molecule(const_mol)
            for dmol in distinct_molecules:
                smiles = stk.Smiles().get_key(dmol)

                if "." in smiles:
                    msg = "Found `.` in smiles."
                    raise RuntimeError(msg)

                if smiles in yielded_smiles:
                    continue
                yielded_smiles.add(smiles)

                # Name based on which bbs are present in how many.
                present_bbs = self.get_present_building_blocks(
                    const_mol=const_mol,
                    subset=dmol,
                )

                intermediate_name = f"idx{len(intermediates)}_"
                for bb, bb_ids in present_bbs.items():
                    num_of_bb = len(bb_ids)
                    num_fgs = bb.get_num_functional_groups()
                    intermediate_name += f"{num_of_bb}x{num_fgs}FG+"
                intermediate_name = intermediate_name[:-1]

                intermediates.append(
                    NamedIntermediate(
                        present_bbs=present_bbs,
                        molecule=dmol,
                        num_atoms=dmol.get_num_atoms(),
                        intermediate_name=intermediate_name,
                        count_bbs=sum(len(i) for i in present_bbs.values()),
                    )
                )

        return intermediates
