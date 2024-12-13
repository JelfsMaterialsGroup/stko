Host-guest complexes with OpenMM
================================

Here we build on the cage optimisation workflow to show how to build host-guest
systems using the open-source `OpenMM <https://openmm.org/>`_ and
`OpenFF <https://openforcefield.org/>`_ infrastructure.

This is shown in an example
`script <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/openmm_hg_example.py>`_
that we run through below.


First we load a cage from the previous example, the classic CC3 porous organic
cage. And we build to `stk.host_guest.Guest` objects.

.. code-block:: python

    from pathlib import Path

    import openmm
    import stk
    from openff.toolkit import ForceField

    import stko

    input_cage = stk.BuildingBlock.init_from_file(
        Path("cage_openmm_example_directory/opt_cage.mol")
    )

    guest1 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("CC"),
        displacement=(0.0, 3.0, 0.0),
    )
    guest2 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("C1CCCC1"),
    )
    hg_complex = stk.ConstructedMolecule(
        topology_graph=stk.host_guest.Complex(
            host=stk.BuildingBlock.init_from_molecule(input_cage),
            guests=(guest1, guest2),
            optimizer=stk.Spinner(),
        ),
    )

.. moldoc::

    import moldoc.molecule as molecule
    import stk

    try:
        hg_complex = stk.BuildingBlock.init_from_file(
            'source/_static/unopt_hgcomplex.mol',
        )
    except OSError:
        hg_complex = stk.BuildingBlock.init_from_file(
            'docs/source/_static/unopt_hgcomplex.mol',
        )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                hg_complex.get_atoms(),
                hg_complex.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in hg_complex.get_bonds()
        ),
    )

Again, we define some settings, this time at lower temperature.

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
    temperature = 300 * openmm.unit.kelvin
    friction = 10 / openmm.unit.picoseconds
    time_step = 1 * openmm.unit.femtoseconds

We can then run an :class:`stk.OptimizerSequence` built from `OpenMM` classes
to get the structure below in a few minutes!

.. code-block:: python

    # Define sequence.
    optimisation_sequence = stko.OptimizerSequence(
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

    optimised_complex = optimisation_sequence.optimize(hg_complex)


.. moldoc::

    import moldoc.molecule as molecule
    import stk

    try:
        hg_complex = stk.BuildingBlock.init_from_file(
            'source/_static/opt_complex.mol',
        )
    except OSError:
        hg_complex = stk.BuildingBlock.init_from_file(
            'docs/source/_static/opt_complex.mol',
        )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                hg_complex.get_atoms(),
                hg_complex.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in hg_complex.get_bonds()
        ),
    )
