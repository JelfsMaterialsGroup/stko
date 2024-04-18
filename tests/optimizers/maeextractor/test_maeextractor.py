import os

from stko import MAEExtractor


def test_maeextractor(case_data):
    """Test :class:`.MAEExtractor`."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    # Extract the lowest energy conformer into its own .mae file.
    extractor = MAEExtractor(
        run_name=os.path.join(dir_path, case_data.run_name),
        n=case_data.n,
    )

    print(extractor.mae_path)
    assert extractor.mae_path.split("/")[-1] == case_data.mae_path

    assert len(extractor.energies) == len(case_data.energies)
    min_energy = 1e24
    for test, known in zip(
        extractor.energies, case_data.energies, strict=True
    ):
        assert test == known
        if test[0] < min_energy:
            min_energy = test[0]

    print(min_energy)
    assert min_energy == extractor.min_energy

    print(extractor.min_energy)
    assert extractor.min_energy == case_data.min_energy

    print(extractor.path)
    assert extractor.path.split("/")[-1] == case_data.path

    os.remove(os.path.join(dir_path, "test-out.mae"))
    if case_data.n == 1:
        os.remove(os.path.join(dir_path, "test-outEXTRACTED_76.mae"))
    elif case_data.n == 4:
        os.remove(os.path.join(dir_path, "test-outEXTRACTED_36_conf_2.mae"))
        os.remove(os.path.join(dir_path, "test-outEXTRACTED_6_conf_3.mae"))
        os.remove(os.path.join(dir_path, "test-outEXTRACTED_76_conf_0.mae"))
        os.remove(os.path.join(dir_path, "test-outEXTRACTED_78_conf_1.mae"))
