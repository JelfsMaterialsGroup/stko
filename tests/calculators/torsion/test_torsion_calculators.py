import pytest
import stk

import stko


def _polymer_angles_match(
    bb1: stk.BuildingBlock,
    bb2: stk.BuildingBlock,
    torsion_index: int = 0,
) -> bool:
    """Attach bb1 to bb2, and test torsion mapping."""
    polymer = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            orientations=(0, 0),
            num_repeating_units=1,
        )
    )

    tors_calculator = stko.MatchedTorsionCalculator()
    tors_results = tors_calculator.get_results(polymer)

    angles = stko.get_torsion_info_angles(
        polymer, list(tors_results.get_torsion_infos())[torsion_index]
    )
    return angles[0] == pytest.approx(angles[1])


def test_torsion_matching() -> None:
    """Confirm torsions are appropriately mapped.

    Check corresponding torsion angles between building blocks and
    constructed molecules to confirm that the angle is the same in
    both.

    """
    bb1 = stk.BuildingBlock("NCCNCCN", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()])
    assert _polymer_angles_match(bb1, bb2)

    # Test case designed such that default rdkit torsion for
    # constructed molecule spans over both building blocks, requiring
    # stko to correctly clean this up.
    bonders3 = (stk.C(1),)
    deleters3 = (stk.H(7),)
    functional_group3 = stk.GenericFunctionalGroup(
        atoms=bonders3 + deleters3, bonders=bonders3, deleters=deleters3
    )
    bb3 = stk.BuildingBlock("CCCC", [functional_group3])
    bonders4 = (stk.C(0),)
    deleters4 = (stk.H(1),)
    functional_group4 = stk.GenericFunctionalGroup(
        atoms=bonders4 + deleters4, bonders=bonders4, deleters=deleters4
    )
    bb4 = stk.BuildingBlock("C", [functional_group4])
    assert _polymer_angles_match(bb3, bb4)

    # Test case for an exterior atom of a building block torsion being
    # deleted in construction.
    bonders5 = (stk.C(1),)
    deleters5 = (stk.C(0), stk.H(4), stk.H(5), stk.H(6))
    functional_group5 = stk.GenericFunctionalGroup(
        atoms=bonders5 + deleters5, bonders=bonders5, deleters=deleters5
    )
    bb5 = stk.BuildingBlock("CCCC", [functional_group5])
    bb6 = bb4
    assert not _polymer_angles_match(bb5, bb6)
