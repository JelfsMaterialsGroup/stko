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
    tors_calculator = stko.TorsionCalculator()
    tors_results = tors_calculator.get_results(bb1)
    torsions = tors_results.get_torsions()
    for t in torsions:
        print(t)

    # Run calculations for constructed molecule.
    tors_calculator = stko.ConstructedMoleculeTorsionCalculator()
    tors_results = tors_calculator.get_results(polymer)
    torsions = tors_results.get_torsions()
    for t in torsions:
        print(t)
    torsion_infos = tors_results.get_torsion_infos()
    for t in torsion_infos:
        print(t)
    print(tors_results.get_atom_maps())


if __name__ == "__main__":
    main()
