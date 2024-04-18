import numpy as np
import pytest
import stko

try:
    import MDAnalysis as mda
except ModuleNotFoundError:
    mda = None


def test_universe(molecule):
    """Test :meth:`.MDAnalysis.get_universe`.

    Parameters
    ----------
    molecule : :class:`stk.Molecule`
        The molecule to test.

    Returns
    -------
    None : :class:`NoneType`

    """
    if mda is None:
        with pytest.raises(stko.WrapperNotInstalledError):
            result = stko.MDAnalysis().get_universe(molecule)
    else:
        result = stko.MDAnalysis().get_universe(molecule)

        assert result.atoms.n_atoms == molecule.get_num_atoms()

        for atom, stk_atom in zip(result.atoms, molecule.get_atoms(), strict=False):
            assert atom.ix == stk_atom.get_id()

        for pos, stk_pos in zip(
            result.atoms.positions, molecule.get_position_matrix(), strict=False
        ):
            assert np.allclose(pos, stk_pos)
