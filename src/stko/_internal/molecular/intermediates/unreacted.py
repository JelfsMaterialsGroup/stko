import itertools as it
import logging
from collections import abc, defaultdict
from dataclasses import dataclass
from functools import partial

import stk
from rdkit import Chem

from stko._internal.molecular.molecular_utilities import separate_molecule

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class NamedIntermediate:
    """Named intermediate container."""

    intermediate_name: str
    molecule: stk.BuildingBlock
    present_bbs: dict[stk.BuildingBlock, list[int]]
    num_atoms: int
    count_bbs: int

    def __str__(self) -> str:
        """String representation."""
        return repr(self)

    def __repr__(self) -> str:
        """String representation."""
        return f"{self.__class__.__name__}({self.intermediate_name})"


@dataclass
class IntermediatePool:
    """Container of a set of intermediates."""

    intermediates: list[NamedIntermediate]

    def __str__(self) -> str:
        """String representation."""
        return repr(self)

    def __repr__(self) -> str:
        """String representation."""
        return f"{self.__class__.__name__}({len(self.intermediates)})"

    def __len__(self) -> int:
        """Get the length based on the number of intermediates."""
        return len(self.intermediates)


class UnreactedTopologyGraph:
    """A class containing topology graphs and performing subset of reactions.

    Use this to get partially reacted topology graphs.

    .. testcode:: unreacted-topology-graph

        import stk
        import stko

        bb1 = stk.BuildingBlock(
            smiles="NCCN", functional_groups=(stk.PrimaryAminoFactory(),)
        )
        bb2 = stk.BuildingBlock(
            smiles="O=CC(C=O)C=O", functional_groups=(stk.AldehydeFactory(),)
        )
        cage_graphs = stko.topology_functions.UnreactedTopologyGraph(
            stk.cage.TwoPlusThree((bb1, bb2))
        )
        # Get a pool of NamedIntermediates with only 1 reaction, which will
        # contain the reacted + the building blocks (there are 2). You can
        # iterate through that pool to get the named intermediate, containing
        # an stk molecule and other information about the intermediate.
        pool = cage_graphs.get_named_intermediates(n=1)

    .. testcode:: unreacted-topology-graph
        :hide:

        assert len(cage_graphs.get_available_reactions()) == 6
        assert len(pool) == 3

    .. moldoc::

        import moldoc.molecule as molecule
        import stk
        import stko

        bb1 = stk.BuildingBlock(
            smiles="NCCN", functional_groups=(stk.PrimaryAminoFactory(),)
        )
        bb2 = stk.BuildingBlock(
            smiles="O=CC(C=O)C=O", functional_groups=(stk.AldehydeFactory(),)
        )
        cage_graphs = stko.topology_functions.UnreactedTopologyGraph(
            stk.cage.TwoPlusThree((bb1, bb2))
        )
        # Get a NamedIntermediate with only 1 reaction, which contains an
        # stk molecule and other information about the intermediate.
        pool = cage_graphs.get_named_intermediates(n=1)

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    pool.intermediates[0].molecule.get_atoms(),
                    pool.intermediates[0].molecule.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in pool.intermediates[0].molecule.get_bonds()
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
    ) -> abc.Iterator[stk.ConstructedMolecule]:
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

    def get_reacted_smiles(self, n: int | None = None) -> set[str]:
        """Yield constructed molecules with n reactions performed."""
        yielded_smiles = set()
        for const_mol in self.yield_constructed_molecules(n=n):
            distinct_molecules = separate_molecule(const_mol)
            for dmol, _ in distinct_molecules:
                smiles = Chem.CanonSmiles(
                    Chem.MolToSmiles(dmol.to_rdkit_mol())
                )

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
        subset_ids: list[int],
    ) -> dict[stk.BuildingBlock, list[int]]:
        """Get the building blocks present in a constructed molecule."""
        bbs: dict[stk.BuildingBlock, list[int]] = defaultdict(list)
        for atom_info in const_mol.get_atom_infos():
            if atom_info.get_atom().get_id() in subset_ids and (
                atom_info.get_building_block_id()
                not in bbs[atom_info.get_building_block()]  # type: ignore[index]
            ):
                bbs[atom_info.get_building_block()].append(  # type: ignore[index]
                    atom_info.get_building_block_id()  # type: ignore[arg-type]
                )

        return bbs

    def get_named_intermediates(
        self, n: int | None = None
    ) -> IntermediatePool:
        """Yield constructed molecules with up to n reactions performed."""
        yielded_smiles = set()
        pool = IntermediatePool(intermediates=[])
        for const_mol in self.yield_constructed_molecules(n=n):
            distinct_molecules = separate_molecule(const_mol)
            for dmol, datom_ids in distinct_molecules:
                smiles = Chem.CanonSmiles(
                    Chem.MolToSmiles(dmol.to_rdkit_mol())
                )

                if "." in smiles:
                    msg = "Found `.` in smiles."
                    raise RuntimeError(msg)

                if smiles in yielded_smiles:
                    continue
                yielded_smiles.add(smiles)

                # Name based on which bbs are present in how many.
                present_bbs = self.get_present_building_blocks(
                    const_mol=const_mol,
                    subset_ids=datom_ids,
                )

                intermediate_name = f"idx{len(pool.intermediates)}_"
                for bb, bb_ids in present_bbs.items():
                    num_of_bb = len(bb_ids)
                    num_fgs = bb.get_num_functional_groups()
                    intermediate_name += f"{num_of_bb}x{num_fgs}FG+"
                intermediate_name = intermediate_name[:-1]

                pool.intermediates.append(
                    NamedIntermediate(
                        present_bbs=present_bbs,
                        molecule=dmol,
                        num_atoms=dmol.get_num_atoms(),
                        intermediate_name=intermediate_name,
                        count_bbs=sum(len(i) for i in present_bbs.values()),
                    )
                )

        return pool
