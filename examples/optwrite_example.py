# ruff: noqa: T201

import logging
from pathlib import Path

import stk

import stko

logger = logging.getLogger(__name__)


def main() -> None:
    """Run the example."""
    bb1 = stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            num_repeating_units=3,
            optimizer=stk.Collapser(scale_steps=False),
        ),
    )

    output_directory = Path("optwrite_output")
    output_directory.mkdir(exist_ok=True, parents=True)

    optimiser = stko.OptWriterSequence(
        optimizers={
            "etkdg": stko.ETKDG(),
            "uff": stko.UFF(),
            "mmff": stko.MMFF(),
        },
        writer=stk.MolWriter(),
        output_directory=output_directory,
    )

    polymer = optimiser.optimize(polymer)

    # Delete some, rerun.
    output_directory.joinpath("mmff_out.mol").unlink()
    polymer = optimiser.optimize(polymer)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
