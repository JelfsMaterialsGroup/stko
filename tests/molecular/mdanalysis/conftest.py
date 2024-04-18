import pytest
import stk

_pd = stk.BuildingBlock(
    smiles="[Pd+2]",
    functional_groups=(stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)),
    position_matrix=[[0, 0, 0]],
)

# Define a bidentate ligand with two functional groups.
_bidentate_ligand = stk.BuildingBlock(
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
_complex = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.CisProtectedSquarePlanar(
        metals=_pd,
        ligands=_bidentate_ligand,
    ),
)


@pytest.fixture(
    params=(
        stk.BuildingBlock("NCCN"),
        stk.BuildingBlock("c1ccccc1"),
        _complex,
    ),
)
def molecule(request):
    """A :class:`.CaseData` instance.
    """
    return request.param
