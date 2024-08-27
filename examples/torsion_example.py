# ruff: noqa: T201
from pathlib import Path

import stk

import stko


def main() -> None:
    """Run the example."""
    bb1 = stk.BuildingBlock("NCCNCCN", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            orientations=[0, 0],
            num_repeating_units=1,
        )
    )

    output_directory = Path("output_directory")
    # Run calculations for bb.
    bb1.write(output_directory / "tors_test_bb1.mol")
    tors_results = stko.TorsionCalculator().get_results(bb1)
    print(tors_results.get_molecule())
    for t, ang in tors_results.get_torsion_angles():
        print(t, ang, t.get_atom_ids())

    # Run calculations for constructed molecule.
    polymer.write(output_directory / "tors_test_polymer.mol")
    tors_results = stko.ConstructedMoleculeTorsionCalculator().get_results(
        polymer
    )
    print(tors_results.get_molecule())
    for torsion, ang in tors_results.get_torsion_angles():
        print(torsion, ang, torsion.get_atom_ids())
    for torsion_info in tors_results.get_torsion_infos():
        print(
            "c",
            torsion_info.get_torsion(),
            torsion_info.get_building_block(),
            torsion_info.get_building_block_id(),
            torsion_info.get_building_block_torsion(),
        )
    print(tors_results.get_torsion_infos_by_building_block())
    for torsion_info in tors_results.get_torsion_infos():
        print(stko.get_torsion_info_angles(polymer, torsion_info))


if __name__ == "__main__":
    main()
