from pathlib import Path

import numpy as np
import pytest
import stk

import stko
from stko import GulpUFFMDOptimizer, GulpUFFOptimizer
from stko._internal.optimizers.utilities import get_metal_atoms


class FakeGulpUFFOptimizer(GulpUFFOptimizer):
    def _check_path(self, path: Path) -> None:  # noqa: ARG002
        return

    def optimize(self, mol: stko.MoleculeT) -> stko.MoleculeT:
        return mol.with_centroid(np.array([1, 3, 3]))


class FakeGulpUFFMDOptimizer(GulpUFFMDOptimizer):
    def _check_path(self, path: Path) -> None:  # noqa: ARG002
        return

    def optimize(self, mol: stko.MoleculeT) -> stko.MoleculeT:
        return mol.with_centroid(np.array([1, 3, 3]))


@pytest.fixture
def position_section() -> str:
    return (
        "\ncartesian\n"
        "C1 core -0.74031 0.03171 0.10194\n"
        "C1 core 0.75974 -0.01804 -0.03434\n"
        "H1 core -1.14779 -0.63142 -0.70939\n"
        "H1 core -1.11274 -0.39237 1.04655\n"
        "H1 core -1.12625 1.0535 -0.10345\n"
        "H1 core 0.97939 0.03507 -1.13219\n"
        "H1 core 1.25573 0.87341 0.40438\n"
        "H1 core 1.13222 -0.95187 0.4265\n"
    )


@pytest.fixture
def bond_section() -> str:
    return (
        "\nconnect 1 2 \n"
        "connect 1 3 \n"
        "connect 1 4 \n"
        "connect 1 5 \n"
        "connect 2 6 \n"
        "connect 2 7 \n"
        "connect 2 8 \n"
    )


@pytest.fixture
def species_section() -> str:
    return "\nspecies\nC1 C_3\nH1 H_\n"


def test_gulp_position_section(
    unoptimized_mol: stk.BuildingBlock,
    position_section: str,
) -> None:
    opt = FakeGulpUFFOptimizer(
        gulp_path="",
        maxcyc=1000,
        metal_FF=None,
        metal_ligand_bond_order=None,
        conjugate_gradient=False,
        output_dir="",
    )
    opt.assign_FF(unoptimized_mol)
    type_translator = opt._type_translator()  # noqa: SLF001
    test = opt._position_section(unoptimized_mol, type_translator)  # noqa: SLF001
    assert position_section == test


def test_gulp_bond_section(
    unoptimized_mol: stk.BuildingBlock,
    bond_section: str,
) -> None:
    opt = FakeGulpUFFOptimizer(gulp_path="")
    opt.assign_FF(unoptimized_mol)
    metal_atoms = get_metal_atoms(unoptimized_mol)
    test = opt._bond_section(unoptimized_mol, metal_atoms)  # noqa: SLF001
    assert bond_section == test


def test_gulp_species_section(
    unoptimized_mol: stk.BuildingBlock,
    species_section: str,
) -> None:
    opt = FakeGulpUFFOptimizer(gulp_path="")
    opt.assign_FF(unoptimized_mol)
    type_translator = opt._type_translator()  # noqa: SLF001
    test = opt._species_section(type_translator)  # noqa: SLF001
    assert species_section == test


@pytest.fixture
def atom_types() -> list[str]:
    return ["C", "C", "H", "H", "H", "H", "H", "H"]


@pytest.fixture
def trajectory() -> dict[int, dict[str, float]]:
    """Defines output of the trajectory properties ignoring coords."""
    return {
        0: {
            "time": 2.99999999999978,
            "KE": 0.148961461675797,
            "E": 0.322914341829668,
            "T": 192.069356724061,
        },
        1: {
            "time": 5.0,
            "KE": 0.207467138253947,
            "E": 0.245970616427967,
            "T": 267.505966560289,
        },
        2: {
            "time": 7.00000000000067,
            "KE": 0.267980085285543,
            "E": 0.207888721719671,
            "T": 345.530730006366,
        },
        3: {
            "time": 9.00000000000045,
            "KE": 0.203650341841703,
            "E": 0.203495426725446,
            "T": 262.584629031782,
        },
        4: {
            "time": 10.9999999999993,
            "KE": 0.209627546107803,
            "E": 0.246948965763409,
            "T": 270.291573938759,
        },
    }


@pytest.fixture
def min_energy_time_step() -> int:
    return 3


@pytest.fixture
def xyz_string() -> str:
    return (
        "8\n"
        "1,5.0,0.207467138253947,0.245970616427967,267.505966560289\n"
        "C -0.74873 0.01254 0.08187\n"
        "C 0.75325 -0.00453 -0.03139\n"
        "H -1.21565 -0.9294 -0.08591\n"
        "H -1.00574 0.41323 1.09345\n"
        "H -1.09608 0.79503 -0.60147\n"
        "H 1.02245 -0.86479 -0.73236\n"
        "H 1.19827 0.93931 -0.43282\n"
        "H 1.25464 -0.29976 0.89505\n"
    )


def test_gulp_convert_traj_to_xyz(
    atom_types: list[str], trajectory: dict[int, dict[str, float]]
) -> None:
    test_dir = Path(__file__).resolve().parent
    test_xyz = Path(f"{test_dir}/fixtures/gulp_MD_template.xyz")
    test_traj = Path(f"{test_dir}/fixtures/gulp_MD.trg")
    opt = FakeGulpUFFMDOptimizer(
        gulp_path="",
    )
    test_atom_types, test_trajectory_data, _ = opt._convert_traj_to_xyz(  # noqa: SLF001
        output_xyz=test_xyz,
        output_traj=test_traj,
    )

    for i, t in zip(atom_types, test_atom_types, strict=False):
        assert i == t

    for ts in trajectory:
        test_ts_dict = test_trajectory_data[ts]
        ts_dict = trajectory[ts]
        assert test_ts_dict["time"] == ts_dict["time"]
        assert test_ts_dict["KE"] == ts_dict["KE"]
        assert test_ts_dict["E"] == ts_dict["E"]
        assert test_ts_dict["T"] == ts_dict["T"]


def test_gulp_calculate_lowest_energy_conformer(
    min_energy_time_step: int,
) -> None:
    test_dir = Path(__file__).resolve().parent
    test_xyz = Path(f"{test_dir}/fixtures/gulp_MD_template.xyz")
    test_traj = Path(f"{test_dir}/fixtures/gulp_MD.trg")
    opt = FakeGulpUFFMDOptimizer(
        gulp_path="",
    )
    atom_types, trajectory_data, xyz_traj_lines = opt._convert_traj_to_xyz(  # noqa: SLF001
        output_xyz=test_xyz,
        output_traj=test_traj,
    )
    min_ts = opt._calculate_lowest_energy_conformer(trajectory_data)  # noqa: SLF001

    assert min_ts == min_energy_time_step


def test_gulp_write_conformer_xyz_file(xyz_string: str) -> None:
    test_dir = Path(__file__).resolve().parent
    test_xyz = Path(f"{test_dir}/fixtures/gulp_MD_template.xyz")
    test_traj = Path(f"{test_dir}/fixtures/gulp_MD.trg")
    opt = FakeGulpUFFMDOptimizer(
        gulp_path="",
    )
    atom_types, trajectory_data, xyz_traj_lines = opt._convert_traj_to_xyz(  # noqa: SLF001
        output_xyz=test_xyz,
        output_traj=test_traj,
    )

    test_output = Path(f"{test_dir}/fixtures/conformer.xyz")
    opt._write_conformer_xyz_file(  # noqa: SLF001
        ts=1,
        ts_data=trajectory_data[1],
        filename=test_output,
        atom_types=atom_types,
    )

    test_string = ""
    with test_output.open() as f:
        for line in f.readlines():
            test_string += line

    assert test_string == xyz_string
