from pytest import approx
import stko
import stk
from examples.torsion_example import get_torsion_info_angles


def test_constructed_molecule():
    def test_polymer(bb1, bb2):
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
        # print(tors_results.get_molecule())
        # print(tors_results.get_torsion_infos_by_building_block())
        for t in tors_results.get_torsion_infos():
            print(f'torsion = {t._torsion}')
            print(f'building_block_torsion = {t._building_block_torsion}')
            print(get_torsion_info_angles(polymer, t))
            break
        angles = get_torsion_info_angles(polymer, list(tors_results.get_torsion_infos())[0])
        assert angles[0] == approx(angles[1])
        
    bb1 = stk.BuildingBlock('NCCNCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    bb3 = stk.BuildingBlock('CCCC', [stk.SingleAtom(stk.C(1))])
    bb4 = stk.BuildingBlock('C', [stk.SingleAtom(stk.C(0))])

    test_polymer(bb1, bb2)
    test_polymer(bb3, bb4)
