from pytest import approx
import stko
import stk


def test_torsion_matching():
    """confirm torsions are appropriately mapped from constructed
    molecules to building blocks
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
        return angles[0] == approx(angles[1])

    bb1 = stk.BuildingBlock('NCCNCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    assert test_polymer(bb1, bb2)

    # test case designed such that default rdkit torsion for
    # constructed molecule spans over both building blocks, requiring
    # stko to correctly clean this up
    bb3 = stk.BuildingBlock('CCCC', [stk.SingleAtom(stk.C(1))])
    bb4 = stk.BuildingBlock('C', [stk.SingleAtom(stk.C(0))])
    assert test_polymer(bb3, bb4)

    # test case for an exterior atom of a building block torsion being
    # deleted in construction
    functional_group = stk.GenericFunctionalGroup(
        atoms=(stk.C(0), stk.C(1)),
        bonders=(stk.C(1),), deleters=(stk.C(0),))
    bb5 = stk.BuildingBlock('CCCC', [functional_group])
    bb6 = stk.BuildingBlock('C', [stk.SingleAtom(stk.C(0))])
    assert not test_polymer(bb5, bb6)
