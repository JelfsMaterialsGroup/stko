from dataclasses import dataclass

import numpy as np
import pytest
import stk

from stko._internal.utilities.exceptions import ExpectedMetalError

from .test_gulp import FakeGulpUFFOptimizer

pd_metal = stk.BuildingBlock(
    smiles="[Pd+2]",
    functional_groups=(stk.SingleAtom(stk.Pd(0, charge=2)) for _ in range(4)),
    position_matrix=np.array([[0, 0, 0]]),
)

# Define a bidentate ligand with two functional groups.
bidentate_ligand = stk.BuildingBlock(
    smiles="NCCN",
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts="[#7]~[#6]",
            bonders=(0,),
            deleters=(),
        ),
    ],
)

# Construct a cis-protected square planar metal complex.
complex_ = stk.ConstructedMolecule(
    stk.metal_complex.BidentateSquarePlanar(
        metals=pd_metal,
        ligands=bidentate_ligand,
    )
)


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    atom_types: dict[int, str]
    has_metal: bool
    metal_ff: dict[int, str] | None


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock("CCC"),
            atom_types={
                0: "C_3",
                1: "C_3",
                2: "C_3",
                3: "H_",
                4: "H_",
                5: "H_",
                6: "H_",
                7: "H_",
                8: "H_",
                9: "H_",
                10: "H_",
            },
            has_metal=False,
            metal_ff=None,
        ),
        CaseData(
            molecule=stk.BuildingBlock.init_from_molecule(complex_),
            atom_types={
                0: "Pd4+2",
                1: "N_3",
                2: "C_3",
                3: "C_3",
                4: "N_3",
                5: "H_",
                6: "H_",
                7: "H_",
                8: "H_",
                9: "H_",
                10: "H_",
                11: "H_",
                12: "H_",
                13: "N_3",
                14: "C_3",
                15: "C_3",
                16: "N_3",
                17: "H_",
                18: "H_",
                19: "H_",
                20: "H_",
                21: "H_",
                22: "H_",
                23: "H_",
                24: "H_",
            },
            has_metal=True,
            metal_ff={46: "Pd4+2"},
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCC"),
            atom_types={
                0: "C_3",
                1: "C_3",
                2: "C_3",
                3: "H_",
                4: "H_",
                5: "H_",
                6: "H_",
                7: "H_",
                8: "H_",
                9: "H_",
                10: "H_",
            },
            has_metal=False,
            metal_ff={26: "Fe4+2"},
        ),
        CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            atom_types={
                0: "C_R",
                1: "C_R",
                2: "C_R",
                3: "C_R",
                4: "C_R",
                5: "C_R",
                6: "H_",
                7: "H_",
                8: "H_",
                9: "H_",
                10: "H_",
                11: "H_",
            },
            has_metal=False,
            metal_ff=None,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                "C1=CC(=CC(=C1)C#CC2=CN=CC=C2)C#CC3=CN=CC=C3"
            ),
            atom_types={
                0: "C_R",
                1: "C_R",
                2: "C_R",
                3: "C_R",
                4: "C_R",
                5: "C_R",
                6: "C_1",
                7: "C_1",
                8: "C_R",
                9: "C_R",
                10: "N_R",
                11: "C_R",
                12: "C_R",
                13: "C_R",
                14: "C_1",
                15: "C_1",
                16: "C_R",
                17: "C_R",
                18: "N_R",
                19: "C_R",
                20: "C_R",
                21: "C_R",
                22: "H_",
                23: "H_",
                24: "H_",
                25: "H_",
                26: "H_",
                27: "H_",
                28: "H_",
                29: "H_",
                30: "H_",
                31: "H_",
                32: "H_",
                33: "H_",
            },
            has_metal=False,
            metal_ff=None,
        ),
        CaseData(
            molecule=stk.BuildingBlock("C1CCC(C(C1)N)N"),
            atom_types={
                0: "C_3",
                1: "C_3",
                2: "C_3",
                3: "C_3",
                4: "C_3",
                5: "C_3",
                6: "N_3",
                7: "N_3",
                8: "H_",
                9: "H_",
                10: "H_",
                11: "H_",
                12: "H_",
                13: "H_",
                14: "H_",
                15: "H_",
                16: "H_",
                17: "H_",
                18: "H_",
                19: "H_",
                20: "H_",
                21: "H_",
            },
            has_metal=False,
            metal_ff=None,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CC=O"),
            atom_types={
                0: "C_3",
                1: "C_2",
                2: "O_2",
                3: "H_",
                4: "H_",
                5: "H_",
                6: "H_",
            },
            has_metal=False,
            metal_ff=None,
        ),
    ),
)
def test_molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param


def test_assign_ff(test_molecule: CaseData) -> None:
    gulp_opt = FakeGulpUFFOptimizer(
        gulp_path="not_required",
        metal_FF=test_molecule.metal_ff,
    )

    # Assign the force field.
    if test_molecule.has_metal and test_molecule.metal_ff is None:
        with pytest.raises(ExpectedMetalError):
            gulp_opt.assign_FF(test_molecule.molecule)
    else:
        gulp_opt.assign_FF(test_molecule.molecule)

        assert len(gulp_opt.atom_labels) == len(test_molecule.atom_types)
        for aid in gulp_opt.atom_labels:
            expected_type = test_molecule.atom_types[aid]
            test = gulp_opt.atom_labels[aid][0]
            assert expected_type == test
