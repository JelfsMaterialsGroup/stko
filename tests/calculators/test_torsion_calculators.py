import stko
import stk


def test_constructed_molecule():
    bb1 = stk.BuildingBlock('CCCC', [stk.SingleAtom(stk.C(1))])
    bb2 = stk.BuildingBlock('C', [stk.SingleAtom(stk.C(0))])

    polymer = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            orientations=[0, 0],
            num_repeating_units=1
        )
    )

    tors_calculator = stko.ConstructedMoleculeTorsionCalculator()
    tors_results = tors_calculator.get_results(polymer)
    # print(tors_results.get_torsion_infos_by_building_block())
    assert True
