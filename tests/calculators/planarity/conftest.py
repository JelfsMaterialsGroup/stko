from dataclasses import dataclass

import numpy as np
import pytest
import stk

import stko

_macrocycle = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
        building_blocks=(
            stk.BuildingBlock(
                smiles="BrCCBr",
                functional_groups=[stk.BromoFactory()],
            ),
            stk.BuildingBlock(
                smiles="BrCNCCBr",
                functional_groups=[stk.BromoFactory()],
            ),
        ),
        repeating_unit="AB",
        num_repeating_units=2,
    ),
)

_square_planar = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.SquarePlanar(
        metals=stk.BuildingBlock(
            smiles="[Pd+2]",
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2)) for _ in range(4)
            ),
            position_matrix=np.array([[0, 0, 0]]),
        ),
        ligands=stk.BuildingBlock(
            smiles="NBr",
            functional_groups=(stk.PrimaryAminoFactory(),),
        ),
        optimizer=stk.MCHammer(),
    )
)
uff = stko.UFF()
_square_planar_uff = uff.optimize(_square_planar)
_octahedral = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.Octahedral(
        metals=stk.BuildingBlock(
            smiles="[Fe+2]",
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2)) for _ in range(6)
            ),
            position_matrix=np.array([[0, 0, 0]]),
        ),
        ligands=stk.BuildingBlock(
            smiles="NBr",
            functional_groups=(stk.PrimaryAminoFactory(),),
        ),
        optimizer=stk.MCHammer(),
    )
)


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    plane_ids: tuple[int, ...] | None
    deviation_ids: tuple[int, ...] | None
    name: str


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            plane_ids=None,
            deviation_ids=None,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)cc3c2c1"
            ),
            plane_ids=None,
            deviation_ids=None,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            plane_ids=None,
            deviation_ids=None,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            plane_ids=None,
            deviation_ids=None,
            name=name,
        ),
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrCCBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles="BrCNCCBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit="AB",
                    num_repeating_units=2,
                ),
            ),
            plane_ids=None,
            deviation_ids=None,
            name=name,
        ),
        lambda name: CaseData(
            molecule=_macrocycle,
            plane_ids=None,
            deviation_ids=None,
            name=name,
        ),
        lambda name: CaseData(
            molecule=_macrocycle,
            plane_ids=tuple(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6  # noqa: PLR2004
            ),
            deviation_ids=tuple(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6  # noqa: PLR2004
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=_macrocycle,
            plane_ids=tuple(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6  # noqa: PLR2004
            ),
            deviation_ids=None,
            name=name,
        ),
        lambda name: CaseData(
            molecule=_macrocycle,
            plane_ids=None,
            deviation_ids=tuple(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6  # noqa: PLR2004
            ),
            name=name,
        ),
        lambda name: CaseData(
            molecule=_square_planar,
            plane_ids=None,
            deviation_ids=None,
            name=name,
        ),
        lambda name: CaseData(
            molecule=_square_planar_uff,
            plane_ids=None,
            deviation_ids=None,
            name=name,
        ),
        lambda name: CaseData(
            molecule=_octahedral,
            plane_ids=None,
            deviation_ids=None,
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
