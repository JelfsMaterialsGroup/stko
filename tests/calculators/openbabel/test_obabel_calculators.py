import pytest
import stk

import stko

try:
    from openbabel import openbabel
except ImportError:
    openbabel = None


def test_open_babel_energy(unoptimized_mol: stk.BuildingBlock) -> None:
    if openbabel is None:
        with pytest.raises(stko.WrapperNotInstalledError):
            calculator = stko.OpenBabelEnergy("uff")
    else:
        calculator = stko.OpenBabelEnergy("uff")
        test_energy = calculator.get_energy(unoptimized_mol)
        assert test_energy == 141.44622279628743  # noqa: PLR2004

        calculator = stko.OpenBabelEnergy("gaff")
        test_energy = calculator.get_energy(unoptimized_mol)
        assert test_energy == 66.47095418890525  # noqa: PLR2004

        calculator = stko.OpenBabelEnergy("ghemical")
        test_energy = calculator.get_energy(unoptimized_mol)
        assert test_energy == 86.59956589041794  # noqa: PLR2004

        calculator = stko.OpenBabelEnergy("mmff94")
        test_energy = calculator.get_energy(unoptimized_mol)
        assert test_energy == 7.607999187460175  # noqa: PLR2004
