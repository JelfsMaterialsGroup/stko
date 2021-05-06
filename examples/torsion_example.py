import stk
import stko


def get_torsion_info_angles(mol, torsion_info):
    torsion = torsion_info.get_torsion()
    angle = stko.calculate_dihedral(
        pt1=tuple(
            mol.get_atomic_positions(
                torsion.get_atom_ids()[0]
            )
        )[0],
        pt2=tuple(
            mol.get_atomic_positions(
                torsion.get_atom_ids()[1]
            )
        )[0],
        pt3=tuple(
            mol.get_atomic_positions(
                torsion.get_atom_ids()[2]
            )
        )[0],
        pt4=tuple(
            mol.get_atomic_positions(
                torsion.get_atom_ids()[3]
            )
        )[0],
    )
    bb_torsion = torsion_info.get_building_block_torsion()
    if bb_torsion is None:
        bb_angle = None
    else:
        bb_angle = stko.calculate_dihedral(
            pt1=tuple(
                torsion_info.get_building_block().get_atomic_positions(
                    bb_torsion.get_atom_ids()[0]
                )
            )[0],
            pt2=tuple(
                torsion_info.get_building_block().get_atomic_positions(
                    bb_torsion.get_atom_ids()[1]
                )
            )[0],
            pt3=tuple(
                torsion_info.get_building_block().get_atomic_positions(
                    bb_torsion.get_atom_ids()[2]
                )
            )[0],
            pt4=tuple(
                torsion_info.get_building_block().get_atomic_positions(
                    bb_torsion.get_atom_ids()[3]
                )
            )[0],
        )
    return angle, bb_angle


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
        print(get_torsion_info_angles(polymer, t))


if __name__ == "__main__":
    main()
