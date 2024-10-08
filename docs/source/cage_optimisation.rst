Cage optimisation workflow
==========================

Here we implement a cage optimisation workflow similar to that found `here <https://pubs.rsc.org/en/content/articlelanding/2018/sc/c8sc03560a>`_
but using the open-source `OpenMM <https://openmm.org/>`_ and
`OpenFF <https://openforcefield.org/>`_ infrastructure.

This is shown in an example
`script <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/cage_openmm_example.py>`_
that we run through below.


First we build a cage, the classic CC3 porous organic cage. But note, we are
not handling the detailed stereochemistry of this system here.

.. code-block:: python

    import openmm
    import stk
    from openff.toolkit import ForceField

    import stko

    # Construct a cage.
    bb1 = stk.BuildingBlock(
        smiles="C1CCC(C(C1)N)N",
        functional_groups=[stk.PrimaryAminoFactory()],
    )
    bb2 = stk.BuildingBlock(
        smiles="C1=C(C=C(C=C1C=O)C=O)C=O",
        functional_groups=[stk.AldehydeFactory()],
    )
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
    )


.. moldoc::

    import moldoc.molecule as molecule
    import stk

    bb1 = stk.BuildingBlock(
        smiles="C1CCC(C(C1)N)N",
        functional_groups=[stk.PrimaryAminoFactory()],
    )
    bb2 = stk.BuildingBlock(
        smiles="C1=C(C=C(C=C1C=O)C=O)C=O",
        functional_groups=[stk.AldehydeFactory()],
    )
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
    )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                cage.get_atoms(),
                cage.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in cage.get_bonds()
        ),
    )

First we define some settings and a function for generating a new integrator
from these settings.

.. code-block:: python

    def integrator(
        *,
        temperature: float,
        friction: float,
        time_step: float,
    ) -> openmm.LangevinIntegrator:
        """Define integrator."""
        integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
        integrator.setRandomNumberSeed(127)
        return integrator

    # Settings.
    force_field = ForceField("openff_unconstrained-2.1.0.offxml")
    partial_charges = "espaloma-am1bcc"
    temperature = 700 * openmm.unit.kelvin
    friction = 10 / openmm.unit.picoseconds
    time_step = 1 * openmm.unit.femtoseconds

We can then run an :class:`stk.OptimizerSequence` built from `OpenMM` classes
to get the structure below in a few minutes!

.. code-block:: python

    # Define sequence.
    optimisation_sequence = stko.OptimizerSequence(
        # Restricted true to optimised the constructed bonds.
        stko.OpenMMForceField(
            force_field=force_field,
            restricted=True,
            partial_charges_method=partial_charges,
        ),
        # Unrestricted optimisation.
        stko.OpenMMForceField(
            # Load the openff-2.1.0 force field appropriate for
            # vacuum calculations (without constraints)
            force_field=force_field,
            restricted=False,
            partial_charges_method=partial_charges,
        ),
        # Molecular dynamics, short for equilibration.
        stko.OpenMMMD(
            force_field=force_field,
            output_directory=output_directory / "md_optimisation",
            integrator=integrator(
                temperature=temperature,
                friction=friction,
                time_step=time_step,
            ),
            random_seed=275,
            partial_charges_method=partial_charges,
            # Frequency here is not related to the num confs tested.
            reporting_freq=100,
            trajectory_freq=100,
            # 10 ps
            num_steps=10_000,
            num_conformers=10,
            platform="CUDA",
            conformer_optimiser=stko.OpenMMForceField(
                force_field=force_field,
                restricted=False,
                partial_charges_method=partial_charges,
            ),
        ),
        # Long MD, for collecting lowest energy conformers.
        stko.OpenMMMD(
            force_field=force_field,
            output_directory=output_directory / "md_optimisation",
            integrator=integrator(
                temperature=temperature,
                friction=friction,
                time_step=time_step,
            ),
            random_seed=275,
            partial_charges_method=partial_charges,
            # Frequency here is not related to the num confs tested.
            reporting_freq=100,
            trajectory_freq=100,
            # 0.2 ns
            num_steps=200_000,
            # 1 every 4 ps
            num_conformers=50,
            platform="CUDA",
            conformer_optimiser=stko.OpenMMForceField(
                force_field=force_field,
                restricted=False,
                partial_charges_method=partial_charges,
            ),
        ),
    )

    optimised_cage = optimisation_sequence.optimize(cage)


.. moldoc::

    import moldoc.molecule as molecule
    import stk
    import numpy as np

    bb1 = stk.BuildingBlock(
        smiles="C1CCC(C(C1)N)N",
        functional_groups=[stk.PrimaryAminoFactory()],
    )
    bb2 = stk.BuildingBlock(
        smiles="C1=C(C=C(C=C1C=O)C=O)C=O",
        functional_groups=[stk.AldehydeFactory()],
    )
    cage = stk.BuildingBlock.init_from_file(
        'source/_static/openmm_opt_file.mol'
    )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                cage.get_atoms(),
                cage.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in cage.get_bonds()
        ),
    )