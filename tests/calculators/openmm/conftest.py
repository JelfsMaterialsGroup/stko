from dataclasses import dataclass

import pytest
import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    energy: float


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            energy=173.9890271340685,
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock("NCCN"),
                    guests=stk.host_guest.Guest(
                        building_block=stk.BuildingBlock("CC"),
                        displacement=(5, 5, 5),
                    ),
                    optimizer=stk.Spinner(),
                ),
            ),
            energy=11783.42209090402,
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock("NCCN"),
                    guests=stk.host_guest.Guest(
                        building_block=stk.BuildingBlock("NCCN"),
                        displacement=(4, 4, 4),
                    ),
                    optimizer=stk.Spinner(),
                ),
            ),
            energy=347.7239171693477,
        ),
    ],
)
def case_molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param
