import stko


def test_open_babel_energy(unoptimized_mol):
    calculator = stko.OpenBabelEnergy('uff')
    test_energy = calculator.get_energy(unoptimized_mol)
    assert test_energy == 141.44622279628743

    calculator = stko.OpenBabelEnergy('gaff')
    test_energy = calculator.get_energy(unoptimized_mol)
    assert test_energy == 66.47095418890525

    calculator = stko.OpenBabelEnergy('ghemical')
    test_energy = calculator.get_energy(unoptimized_mol)
    assert test_energy == 86.59956589041794

    calculator = stko.OpenBabelEnergy('mmff94')
    test_energy = calculator.get_energy(unoptimized_mol)
    assert test_energy == 7.607999187460175
