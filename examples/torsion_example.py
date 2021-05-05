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
    for t, ang in zip(
        tors_results.get_torsions(),
        tors_results.get_torsion_angles()
    ):
        print(t, ang, t.get_atom_ids())

    # Run calculations for constructed molecule.
    polymer.write('output_directory/tors_test_polymer.mol')
    tors_calculator = stko.ConstructedMoleculeTorsionCalculator()
    tors_results = tors_calculator.get_results(polymer)
    print(tors_results.get_molecule())
    for t, ang in zip(
        tors_results.get_torsions(),
        tors_results.get_torsion_angles()
    ):
        print(t, ang, t.get_atom_ids())
    torsion_infos = tors_results.get_torsion_infos()
    for t in torsion_infos:
        print(t)
    print(tors_results.get_atom_maps())


if __name__ == "__main__":
    main()
