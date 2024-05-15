import numpy as np
import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            constructed_molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()]),
                        stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()]),
                    ),
                    repeating_unit="AB",
                    num_repeating_units=2,
                ),
            ),
            centroids={
                0: np.array([-0.27896625, 0.78273808, -0.04626922]),
                1: np.array([12.7437916, 1.08096657, 0.04936503]),
                2: np.array([6.39811868, -0.2408838, -0.83794888]),
                3: np.array([19.24738842, -0.13160962, -0.70954977]),
            },
            atom_ids={
                0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                1: [10, 11, 12, 13, 14, 15, 16, 17],
                2: [18, 19, 20, 21, 22, 23, 24, 25, 26, 27],
                3: [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38],
            },
            name=name,
        ),
        lambda name: CaseData(
            constructed_molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()]),
                        stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()]),
                    ),
                    repeating_unit="AB",
                    num_repeating_units=2,
                    optimizer=stk.Collapser(scale_steps=False),
                ),
            ),
            centroids={
                0: np.array([3.12249409, 0.62482436, -0.17366926]),
                1: np.array([11.66672552, 0.82049207, -0.11092363]),
                2: np.array([7.50332952, -0.04677396, -0.69309028]),
                3: np.array([15.93373539, 0.02492083, -0.60884763]),
            },
            atom_ids={
                0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                1: [10, 11, 12, 13, 14, 15, 16, 17],
                2: [18, 19, 20, 21, 22, 23, 24, 25, 26, 27],
                3: [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38],
            },
            name=name,
        ),
        lambda name: CaseData(
            constructed_molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    (
                        stk.BuildingBlock(
                            smiles="NCCN",
                            functional_groups=[stk.PrimaryAminoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles="O=CC(C=O)C=O",
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                    )
                ),
            ),
            centroids={
                0: np.array([0.0952583, 0.05782588, 7.70819073]),
                1: np.array([-6.26032559, -3.61722913, -2.67443437]),
                2: np.array([6.26277511, -3.61298643, -2.67443437]),
                3: np.array([-2.44952914e-03, 7.23021556e00, -2.67443437e00]),
                4: np.array([-4.70844968, -2.71842469, 3.85903949]),
                5: np.array([4.71605668, -2.72476131, 3.85903947]),
                6: np.array([-0.00924248, 5.44294872, 3.85903889]),
                7: np.array([-0.01403689, -5.46061179, -3.84785617]),
                8: np.array([-4.73393618, 2.72180559, -3.84785617]),
                9: np.array([4.73604697, 2.71814959, -3.84785617]),
            },
            atom_ids={
                0: [0, 1, 2, 3, 4, 5, 6, 7],
                1: [8, 9, 10, 11, 12, 13, 14, 15],
                2: [16, 17, 18, 19, 20, 21, 22, 23],
                3: [24, 25, 26, 27, 28, 29, 30, 31],
                4: [32, 33, 34, 35, 36, 37, 38, 39],
                5: [40, 41, 42, 43, 44, 45, 46, 47],
                6: [48, 49, 50, 51, 52, 53, 54, 55],
                7: [56, 57, 58, 59, 60, 61, 62, 63],
                8: [64, 65, 66, 67, 68, 69, 70, 71],
                9: [72, 73, 74, 75, 76, 77, 78, 79],
            },
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
