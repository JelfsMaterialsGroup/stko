import numpy as np
import pytest
import stk

import stko

try:
    import MDAnalysis as mda  # noqa: N813
except ModuleNotFoundError:
    mda = None


def test_universe(molecule: stk.Molecule) -> None:
    if mda is None:
        with pytest.raises(stko.WrapperNotInstalledError):
            result = stko.MDAnalysis().get_universe(molecule)
    else:
        result = stko.MDAnalysis().get_universe(molecule)

        assert result.atoms.n_atoms == molecule.get_num_atoms()

        for atom, stk_atom in zip(
            result.atoms, molecule.get_atoms(), strict=False
        ):
            assert atom.ix == stk_atom.get_id()

        for pos, stk_pos in zip(
            result.atoms.positions,
            molecule.get_position_matrix(),
            strict=False,
        ):
            assert np.allclose(pos, stk_pos)
