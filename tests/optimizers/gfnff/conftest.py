import numpy as np
import pytest
import stk


@pytest.fixture
def polymer() -> stk.ConstructedMolecule:
    building_block = stk.BuildingBlock(
        smiles="BrCCBr",
        functional_groups=[stk.BromoFactory()],
    )
    return stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(building_block,),
            repeating_unit="A",
            num_repeating_units=3,
        ),
    )


@pytest.fixture
def lantern() -> stk.ConstructedMolecule:
    palladium = stk.BuildingBlock(
        smiles="[Pd+2]",
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2)) for _ in range(4)
        ),
        position_matrix=np.array([[0.0, 0.0, 0.0]]),
    )
    ligand = stk.BuildingBlock(
        smiles="C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )
    return stk.ConstructedMolecule(
        topology_graph=stk.cage.M2L4Lantern(
            building_blocks=(palladium, ligand),
        ),
    )
