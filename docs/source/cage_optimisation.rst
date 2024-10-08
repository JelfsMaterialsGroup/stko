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
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
    )

    # Load this to avoid rerunning!
    new_position_matrix = np.array(
        [
            [-1.97920e00, -2.33840e00, 5.29330e00],
            [-3.21850e00, -2.01280e00, 4.70720e00],
            [-3.34180e00, -8.08100e-01, 3.98840e00],
            [-2.24590e00, 6.33000e-02, 3.84380e00],
            [-1.01290e00, -2.78700e-01, 4.43400e00],
            [-8.69100e-01, -1.47970e00, 5.15870e00],
            [4.53900e-01, -1.83740e00, 5.72770e00],
            [-2.40660e00, 1.29330e00, 3.03360e00],
            [-4.35490e00, -2.95890e00, 4.81740e00],
            [-1.84870e00, -3.26990e00, 5.82810e00],
            [-4.28960e00, -5.82700e-01, 3.51760e00],
            [-1.72400e-01, 3.87100e-01, 4.29010e00],
            [1.27980e00, -1.16640e00, 5.53820e00],
            [-3.36960e00, 1.46290e00, 2.57330e00],
            [-4.18980e00, -3.89170e00, 5.33820e00],
            [-4.78030e00, -5.20410e00, -1.51600e-01],
            [-4.32470e00, -6.05420e00, 8.75300e-01],
            [-3.35590e00, -7.03070e00, 5.80200e-01],
            [-2.83840e00, -7.16300e00, -7.22200e-01],
            [-3.32370e00, -6.32430e00, -1.74500e00],
            [-4.29220e00, -5.34060e00, -1.46720e00],
            [-4.76210e00, -4.41300e00, -2.52130e00],
            [-1.74530e00, -8.11840e00, -1.00840e00],
            [-4.81250e00, -5.91990e00, 2.26630e00],
            [-5.50320e00, -4.43770e00, 9.56000e-02],
            [-2.96510e00, -7.67180e00, 1.35900e00],
            [-2.93830e00, -6.39680e00, -2.75290e00],
            [-5.12880e00, -3.45020e00, -2.19650e00],
            [-1.15780e00, -7.93860e00, -1.89730e00],
            [-4.31490e00, -6.49550e00, 3.03550e00],
            [2.34700e00, -5.62450e00, 2.47440e00],
            [3.02530e00, -4.39100e00, 2.40970e00],
            [3.47730e00, -3.91840e00, 1.16160e00],
            [3.23420e00, -4.65320e00, -1.46000e-02],
            [2.55810e00, -5.88630e00, 6.89000e-02],
            [2.11900e00, -6.38440e00, 1.31040e00],
            [1.39210e00, -7.67220e00, 1.41080e00],
            [3.62010e00, -4.12420e00, -1.34110e00],
            [3.22140e00, -3.56040e00, 3.62070e00],
            [1.98180e00, -5.96360e00, 3.43530e00],
            [3.98010e00, -2.96410e00, 1.07840e00],
            [2.34760e00, -6.46220e00, -8.21700e-01],
            [8.72800e-01, -7.89220e00, 2.33430e00],
            [3.06340e00, -4.48890e00, -2.19170e00],
            [3.59440e00, -2.55290e00, 3.49350e00],
            [3.00000e-04, 1.77600e-01, -1.65950e00],
            [-1.40060e00, 2.90700e-01, -1.74110e00],
            [-2.10490e00, -5.31900e-01, -2.64140e00],
            [-1.42250e00, -1.45750e00, -3.45560e00],
            [-1.76000e-02, -1.54710e00, -3.37100e00],
            [7.00700e-01, -7.33600e-01, -2.47340e00],
            [2.17410e00, -8.22100e-01, -2.34730e00],
            [-2.15950e00, -2.36940e00, -4.36030e00],
            [-2.14360e00, 1.23390e00, -8.73600e-01],
            [5.26300e-01, 7.92300e-01, -9.40000e-01],
            [-3.18380e00, -4.86800e-01, -2.70140e00],
            [5.30900e-01, -2.25820e00, -3.97360e00],
            [2.63710e00, -3.12700e-01, -1.51260e00],
            [-1.62130e00, -3.20460e00, -4.78840e00],
            [-3.22280e00, 1.16310e00, -8.48300e-01],
            [-8.74750e00, -5.40980e00, 3.24400e00],
            [-9.06590e00, -3.95330e00, 3.63770e00],
            [-7.80370e00, -3.06950e00, 3.54650e00],
            [-6.66100e00, -3.63230e00, 4.40950e00],
            [-6.31390e00, -5.09780e00, 3.97920e00],
            [-7.57910e00, -5.97870e00, 4.08540e00],
            [-5.84770e00, -5.16980e00, 2.57580e00],
            [-5.53620e00, -2.67670e00, 4.30830e00],
            [-8.47810e00, -5.44660e00, 2.18510e00],
            [-9.63400e00, -6.03340e00, 3.36800e00],
            [-9.45770e00, -3.92450e00, 4.65600e00],
            [-9.84120e00, -3.55490e00, 2.98070e00],
            [-8.03370e00, -2.05360e00, 3.86880e00],
            [-7.46610e00, -3.00740e00, 2.51010e00],
            [-6.99990e00, -3.64610e00, 5.44680e00],
            [-5.55320e00, -5.49450e00, 4.65370e00],
            [-7.34530e00, -6.98880e00, 3.74590e00],
            [-7.88130e00, -6.05350e00, 5.13020e00],
            [4.35010e00, -5.02540e00, 7.27790e00],
            [3.14570e00, -5.17030e00, 8.23060e00],
            [1.83780e00, -4.71390e00, 7.54810e00],
            [1.95040e00, -3.27500e00, 7.01310e00],
            [3.14870e00, -3.15070e00, 6.01020e00],
            [4.45410e00, -3.58880e00, 6.71120e00],
            [2.96360e00, -4.02000e00, 4.82660e00],
            [6.24200e-01, -2.92360e00, 6.45310e00],
            [4.24210e00, -5.72990e00, 6.44920e00],
            [5.27180e00, -5.28810e00, 7.79930e00],
            [3.31460e00, -4.57440e00, 9.12910e00],
            [3.04910e00, -6.21070e00, 8.54680e00],
            [1.00840e00, -4.77320e00, 8.25350e00],
            [1.59720e00, -5.38120e00, 6.71820e00],
            [2.13990e00, -2.61690e00, 7.86240e00],
            [3.25030e00, -2.10820e00, 5.70470e00],
            [5.28010e00, -3.53660e00, 6.00040e00],
            [4.68100e00, -2.89180e00, 7.51830e00],
            [-1.02780e00, 5.10990e00, -3.06000e-01],
            [-3.33300e-01, 5.33370e00, 1.05260e00],
            [-2.36700e-01, 4.01420e00, 1.84980e00],
            [-1.62180e00, 3.37150e00, 2.03940e00],
            [-2.30760e00, 3.11200e00, 6.52600e-01],
            [-2.41170e00, 4.44160e00, -1.25400e-01],
            [-1.52870e00, 2.17050e00, -1.83300e-01],
            [-1.43050e00, 2.16190e00, 2.87270e00],
            [-3.98700e-01, 4.47130e00, -9.31800e-01],
            [-1.13890e00, 6.06020e00, -8.30400e-01],
            [-8.92800e-01, 6.06900e00, 1.63340e00],
            [6.67000e-01, 5.74010e00, 8.92400e-01],
            [2.12600e-01, 4.19820e00, 2.82620e00],
            [4.12900e-01, 3.30900e00, 1.32780e00],
            [-2.23900e00, 4.07250e00, 2.60390e00],
            [-3.31510e00, 2.72810e00, 8.18700e-01],
            [-2.85780e00, 4.25330e00, -1.10310e00],
            [-3.07840e00, 5.11980e00, 4.07700e-01],
            [2.38150e00, -1.10436e01, -9.31800e-01],
            [1.28440e00, -1.11927e01, -2.00940e00],
            [2.18500e-01, -1.00819e01, -1.87540e00],
            [-3.83200e-01, -1.00826e01, -4.47300e-01],
            [7.22600e-01, -9.86240e00, 6.20900e-01],
            [1.77920e00, -1.09799e01, 4.94200e-01],
            [1.42250e00, -8.56620e00, 4.44300e-01],
            [-1.50840e00, -9.14190e00, -2.11800e-01],
            [2.94940e00, -1.01287e01, -1.11960e00],
            [3.08630e00, -1.18738e01, -9.99000e-01],
            [8.06500e-01, -1.21694e01, -1.91270e00],
            [1.73630e00, -1.11510e01, -3.00230e00],
            [-5.72600e-01, -1.02374e01, -2.60990e00],
            [6.82000e-01, -9.11500e00, -2.07910e00],
            [-7.92000e-01, -1.10790e01, -2.72400e-01],
            [2.64500e-01, -9.92050e00, 1.61100e00],
            [2.57430e00, -1.08053e01, 1.22040e00],
            [1.32150e00, -1.19378e01, 7.42900e-01],
            [-5.76750e00, -1.34190e00, -6.35280e00],
            [-6.91450e00, -2.04920e00, -5.59810e00],
            [-6.37120e00, -2.84330e00, -4.39030e00],
            [-5.29020e00, -3.85340e00, -4.85000e00],
            [-4.11410e00, -3.13440e00, -5.56380e00],
            [-4.65930e00, -2.34010e00, -6.77120e00],
            [-3.42630e00, -2.17840e00, -4.66430e00],
            [-4.79070e00, -4.77070e00, -3.79100e00],
            [-5.33080e00, -5.74500e-01, -5.70860e00],
            [-6.16030e00, -8.32600e-01, -7.23420e00],
            [-7.43810e00, -2.72640e00, -6.27540e00],
            [-7.64070e00, -1.30960e00, -5.25590e00],
            [-7.18390e00, -3.37660e00, -3.89570e00],
            [-5.94310e00, -2.14310e00, -3.67170e00],
            [-5.75710e00, -4.49780e00, -5.59690e00],
            [-3.41620e00, -3.89280e00, -5.92580e00],
            [-3.83900e00, -1.80050e00, -7.24690e00],
            [-5.05640e00, -3.03680e00, -7.50990e00],
            [4.77970e00, -1.76860e00, -5.59980e00],
            [5.36040e00, -3.16920e00, -5.30700e00],
            [4.76900e00, -3.74410e00, -4.00200e00],
            [5.02460e00, -2.77390e00, -2.82240e00],
            [4.40440e00, -1.37650e00, -3.09360e00],
            [4.98680e00, -8.05200e-01, -4.40520e00],
            [2.92930e00, -1.41850e00, -3.24700e00],
            [4.62930e00, -3.28780e00, -1.48570e00],
            [3.70980e00, -1.85690e00, -5.80390e00],
            [5.24190e00, -1.35220e00, -6.49640e00],
            [6.44700e00, -3.10690e00, -5.22070e00],
            [5.13930e00, -3.84090e00, -6.13850e00],
            [5.21660e00, -4.71420e00, -3.78280e00],
            [3.69660e00, -3.89140e00, -4.13610e00],
            [6.10420e00, -2.62450e00, -2.76710e00],
            [4.68230e00, -7.12500e-01, -2.27190e00],
            [4.51020e00, 1.52300e-01, -4.62100e00],
            [6.05170e00, -6.12100e-01, -4.27230e00],
        ]
    )
    cage = cage.with_position_matrix(new_position_matrix)

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