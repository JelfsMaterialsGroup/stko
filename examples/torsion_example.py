import stk
import stko


def main():
    bb1 = stk.BuildingBlock('NCCNCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            orientations=[0, 0],
            num_repeating_units=1
        )
    )

    # Run calculations for bb.
    bb1.write('output_directory/tors_test_bb1.mol')
    tors_calculator = stko.TorsionCalculator()
    tors_results = tors_calculator.get_results(bb1)
    print(tors_results.get_molecule())
    for t, ang in tors_results.get_torsion_angles():
        print(t, ang, t.get_atom_ids())

    # Run calculations for constructed molecule.
    polymer.write('output_directory/tors_test_polymer.mol')
    tors_calculator = stko.ConstructedMoleculeTorsionCalculator()
    tors_results = tors_calculator.get_results(polymer)
    print(tors_results.get_molecule())
    for t, ang in tors_results.get_torsion_angles():
        print(t, ang, t.get_atom_ids())
    for t in tors_results.get_torsion_infos():
        print(
            'c', t.get_torsion(),
            t.get_building_block(),
            t.get_building_block_id(),
            t.get_building_block_torsion(),
        )
    print(tors_results.get_torsion_infos_by_building_block())
    for t in tors_results.get_torsion_infos():
        print(stko.get_torsion_info_angles(polymer, t))


if __name__ == "__main__":
    main()
