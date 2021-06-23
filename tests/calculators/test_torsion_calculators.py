from pytest import approx
import stko
import stk


def test_constructed_molecule():
    """confirm torsions are appropriately mapped from constructed molecules to
    building blocks
    """

    def test_polymer(bb1: stk.BuildingBlock, bb2: stk.BuildingBlock):
        """attach bb1 to bb2 and test torsion mapping"""
        polymer = stk.ConstructedMolecule(
            stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit="AB",
                orientations=[0, 0],
                num_repeating_units=1
            )
        )

        tors_calculator = stko.MatchedTorsionCalculator()
        tors_results = tors_calculator.get_results(polymer)

        angles = stko.get_torsion_info_angles(
            polymer, list(tors_results.get_torsion_infos())[0])
        assert angles[0] == approx(angles[1])

    bb1 = stk.BuildingBlock('NCCNCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    test_polymer(bb1, bb2)

    # test case designed such that default rdkit torsion for constructed
    # molecule spans over both building blocks, requiring stko to
    # correctly clean this up
    bb3 = stk.BuildingBlock('CCCC', [stk.SingleAtom(stk.C(1))])
    bb4 = stk.BuildingBlock('C', [stk.SingleAtom(stk.C(0))])
    test_polymer(bb3, bb4)
