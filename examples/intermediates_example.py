import logging
from pathlib import Path

import stk

import stko

logger = logging.getLogger(__name__)


def main() -> None:
    """Run the example."""
    # Two building blocks you want to react to form a topology graph.
    bb1 = stk.BuildingBlock(
        smiles="NCCN", functional_groups=(stk.PrimaryAminoFactory(),)
    )
    bb2 = stk.BuildingBlock(
        smiles="O=CC(C=O)C=O", functional_groups=(stk.AldehydeFactory(),)
    )

    examples_output = Path("out_intermediates")
    examples_output.mkdir(parents=True, exist_ok=True)

    # Produce a `TopologyGraph` without doing any reactions.
    cage_graphs = stko.topology_functions.UnreactedTopologyGraph(
        stk.cage.FourPlusSix((bb1, bb2))
    )
    logger.info(
        "there are %s possible reactions",
        len(cage_graphs.get_available_reactions()),
    )

    # With up to N reactions performed.
    intermediate_pool = cage_graphs.get_named_intermediates(n=4)
    logger.info("there are %s structures with n=%s", len(intermediate_pool), 4)
    for named_intermediate in intermediate_pool.intermediates:
        named_intermediate.molecule.write(
            examples_output / f"{named_intermediate.intermediate_name}.mol"
        )

    # Now iterate over all possible reactions with varying amounts of
    # completeness and get their smiles.
    all_possible_smiles = cage_graphs.get_reacted_smiles()
    logger.info("there are %s unique smiles", len(all_possible_smiles))


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
