try:
    from openff.toolkit import ForceField
except ImportError:
    ForceField = None
import numpy as np

import stko

from .conftest import CaseData


def test_openmm(case_molecule: CaseData) -> None:
    if ForceField is not None:
        optimiser = stko.OpenMMForceField(
            # Load the openff-2.1.0 force field appropriate for
            # vacuum calculations (without constraints)
            force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
            restricted=False,
            partial_charges_method="espaloma-am1bcc",
        )
        opt_molecule = optimiser.optimize(case_molecule.molecule)

        energy = (
            stko.OpenMMEnergy(
                force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
                partial_charges_method="espaloma-am1bcc",
            ).get_energy(opt_molecule),
        )

        assert np.isclose(energy, case_molecule.opt_energy)
