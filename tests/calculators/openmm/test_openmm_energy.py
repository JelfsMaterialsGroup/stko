try:
    from openff.toolkit import ForceField
except ImportError:
    ForceField = None
import numpy as np

import stko

from .conftest import CaseData


def test_openmm_energy(case_molecule: CaseData) -> None:
    if ForceField is not None:
        energy = (
            stko.OpenMMEnergy(
                force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
                partial_charges_method="espaloma-am1bcc",
            ).get_energy(case_molecule.molecule),
        )

        assert np.isclose(energy, case_molecule.energy)
