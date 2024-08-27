from pathlib import Path

from stko import MAEExtractor
from tests.optimizers.maeextractor.case_data import CaseData


def test_maeextractor(case_data: CaseData) -> None:
    """Test :class:`.MAEExtractor`."""
    dir_path = Path(__file__).resolve().parent
    # Extract the lowest energy conformer into its own .mae file.
    extractor = MAEExtractor(
        run_name=str(dir_path / case_data.run_name),
        n=case_data.n,
    )

    assert extractor.mae_path.name == case_data.mae_path

    assert len(extractor.energies) == len(case_data.energies)
    min_energy = 1e24
    for test, known in zip(
        extractor.energies, case_data.energies, strict=True
    ):
        assert test == known
        min_energy = min(test[0], min_energy)

    assert min_energy == extractor.min_energy

    assert extractor.min_energy == case_data.min_energy

    assert extractor.path.name == case_data.path

    dir_path.joinpath("test-out.mae").unlink()
    if case_data.n == 1:
        dir_path.joinpath("test-outEXTRACTED_76.mae").unlink()
    elif case_data.n == 4:  # noqa: PLR2004
        dir_path.joinpath("test-outEXTRACTED_36_conf_2.mae").unlink()
        dir_path.joinpath("test-outEXTRACTED_6_conf_3.mae").unlink()
        dir_path.joinpath("test-outEXTRACTED_76_conf_0.mae").unlink()
        dir_path.joinpath("test-outEXTRACTED_78_conf_1.mae").unlink()
