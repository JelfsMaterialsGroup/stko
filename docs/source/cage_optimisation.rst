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
    cage = cage.with_position_matrix(
        np.array(
            (
                (-1.9792, -2.3384, 5.2933),
                (-3.2185, -2.0128, 4.7072),
                (-3.3418, -0.8081, 3.9884),
                (-2.2459, 0.0633, 3.8438),
                (-1.0129, -0.2787, 4.4340),
                (-0.8691, -1.4797, 5.1587),
                (0.4539, -1.8374, 5.7277),
                (-2.4066, 1.2933, 3.0336),
                (-4.3549, -2.9589, 4.8174),
                (-1.8487, -3.2699, 5.8281),
                (-4.2896, -0.5827, 3.5176),
                (-0.1724, 0.3871, 4.2901),
                (1.2798, -1.1664, 5.5382),
                (-3.3696, 1.4629, 2.5733),
                (-4.1898, -3.8917, 5.3382),
                (-4.7803, -5.2041, -0.1516),
                (-4.3247, -6.0542, 0.8753),
                (-3.3559, -7.0307, 0.5802),
                (-2.8384, -7.1630, -0.7222),
                (-3.3237, -6.3243, -1.7450),
                (-4.2922, -5.3406, -1.4672),
                (-4.7621, -4.4130, -2.5213),
                (-1.7453, -8.1184, -1.0084),
                (-4.8125, -5.9199, 2.2663),
                (-5.5032, -4.4377, 0.0956),
                (-2.9651, -7.6718, 1.3590),
                (-2.9383, -6.3968, -2.7529),
                (-5.1288, -3.4502, -2.1965),
                (-1.1578, -7.9386, -1.8973),
                (-4.3149, -6.4955, 3.0355),
                (2.3470, -5.6245, 2.4744),
                (3.0253, -4.3910, 2.4097),
                (3.4773, -3.9184, 1.1616),
                (3.2342, -4.6532, -0.0146),
                (2.5581, -5.8863, 0.0689),
                (2.1190, -6.3844, 1.3104),
                (1.3921, -7.6722, 1.4108),
                (3.6201, -4.1242, -1.3411),
                (3.2214, -3.5604, 3.6207),
                (1.9818, -5.9636, 3.4353),
                (3.9801, -2.9641, 1.0784),
                (2.3476, -6.4622, -0.8217),
                (0.8728, -7.8922, 2.3343),
                (3.0634, -4.4889, -2.1917),
                (3.5944, -2.5529, 3.4935),
                (0.0003, 0.1776, -1.6595),
                (-1.4006, 0.2907, -1.7411),
                (-2.1049, -0.5319, -2.6414),
                (-1.4225, -1.4575, -3.4556),
                (-0.0176, -1.5471, -3.3710),
                (0.7007, -0.7336, -2.4734),
                (2.1741, -0.8221, -2.3473),
                (-2.1595, -2.3694, -4.3603),
                (-2.1436, 1.2339, -0.8736),
                (0.5263, 0.7923, -0.9400),
                (-3.1838, -0.4868, -2.7014),
                (0.5309, -2.2582, -3.9736),
                (2.6371, -0.3127, -1.5126),
                (-1.6213, -3.2046, -4.7884),
                (-3.2228, 1.1631, -0.8483),
                (-8.7475, -5.4098, 3.2440),
                (-9.0659, -3.9533, 3.6377),
                (-7.8037, -3.0695, 3.5465),
                (-6.6610, -3.6323, 4.4095),
                (-6.3139, -5.0978, 3.9792),
                (-7.5791, -5.9787, 4.0854),
                (-5.8477, -5.1698, 2.5758),
                (-5.5362, -2.6767, 4.3083),
                (-8.4781, -5.4466, 2.1851),
                (-9.6340, -6.0334, 3.3680),
                (-9.4577, -3.9245, 4.6560),
                (-9.8412, -3.5549, 2.9807),
                (-8.0337, -2.0536, 3.8688),
                (-7.4661, -3.0074, 2.5101),
                (-6.9999, -3.6461, 5.4468),
                (-5.5532, -5.4945, 4.6537),
                (-7.3453, -6.9888, 3.7459),
                (-7.8813, -6.0535, 5.1302),
                (4.3501, -5.0254, 7.2779),
                (3.1457, -5.1703, 8.2306),
                (1.8378, -4.7139, 7.5481),
                (1.9504, -3.2750, 7.0131),
                (3.1487, -3.1507, 6.0102),
                (4.4541, -3.5888, 6.7112),
                (2.9636, -4.0200, 4.8266),
                (0.6242, -2.9236, 6.4531),
                (4.2421, -5.7299, 6.4492),
                (5.2718, -5.2881, 7.7993),
                (3.3146, -4.5744, 9.1291),
                (3.0491, -6.2107, 8.5468),
                (1.0084, -4.7732, 8.2535),
                (1.5972, -5.3812, 6.7182),
                (2.1399, -2.6169, 7.8624),
                (3.2503, -2.1082, 5.7047),
                (5.2801, -3.5366, 6.0004),
                (4.6810, -2.8918, 7.5183),
                (-1.0278, 5.1099, -0.3060),
                (-0.3333, 5.3337, 1.0526),
                (-0.2367, 4.0142, 1.8498),
                (-1.6218, 3.3715, 2.0394),
                (-2.3076, 3.1120, 0.6526),
                (-2.4117, 4.4416, -0.1254),
                (-1.5287, 2.1705, -0.1833),
                (-1.4305, 2.1619, 2.8727),
                (-0.3987, 4.4713, -0.9318),
                (-1.1389, 6.0602, -0.8304),
                (-0.8928, 6.0690, 1.6334),
                (0.6670, 5.7401, 0.8924),
                (0.2126, 4.1982, 2.8262),
                (0.4129, 3.3090, 1.3278),
                (-2.2390, 4.0725, 2.6039),
                (-3.3151, 2.7281, 0.8187),
                (-2.8578, 4.2533, -1.1031),
                (-3.0784, 5.1198, 0.4077),
                (2.3815, -11.0436, -0.9318),
                (1.2844, -11.1927, -2.0094),
                (0.2185, -10.0819, -1.8754),
                (-0.3832, -10.0826, -0.4473),
                (0.7226, -9.8624, 0.6209),
                (1.7792, -10.9799, 0.4942),
                (1.4225, -8.5662, 0.4443),
                (-1.5084, -9.1419, -0.2118),
                (2.9494, -10.1287, -1.1196),
                (3.0863, -11.8738, -0.9990),
                (0.8065, -12.1694, -1.9127),
                (1.7363, -11.1510, -3.0023),
                (-0.5726, -10.2374, -2.6099),
                (0.6820, -9.1150, -2.0791),
                (-0.7920, -11.0790, -0.2724),
                (0.2645, -9.9205, 1.6110),
                (2.5743, -10.8053, 1.2204),
                (1.3215, -11.9378, 0.7429),
                (-5.7675, -1.3419, -6.3528),
                (-6.9145, -2.0492, -5.5981),
                (-6.3712, -2.8433, -4.3903),
                (-5.2902, -3.8534, -4.8500),
                (-4.1141, -3.1344, -5.5638),
                (-4.6593, -2.3401, -6.7712),
                (-3.4263, -2.1784, -4.6643),
                (-4.7907, -4.7707, -3.7910),
                (-5.3308, -0.5745, -5.7086),
                (-6.1603, -0.8326, -7.2342),
                (-7.4381, -2.7264, -6.2754),
                (-7.6407, -1.3096, -5.2559),
                (-7.1839, -3.3766, -3.8957),
                (-5.9431, -2.1431, -3.6717),
                (-5.7571, -4.4978, -5.5969),
                (-3.4162, -3.8928, -5.9258),
                (-3.8390, -1.8005, -7.2469),
                (-5.0564, -3.0368, -7.5099),
                (4.7797, -1.7686, -5.5998),
                (5.3604, -3.1692, -5.3070),
                (4.7690, -3.7441, -4.0020),
                (5.0246, -2.7739, -2.8224),
                (4.4044, -1.3765, -3.0936),
                (4.9868, -0.8052, -4.4052),
                (2.9293, -1.4185, -3.2470),
                (4.6293, -3.2878, -1.4857),
                (3.7098, -1.8569, -5.8039),
                (5.2419, -1.3522, -6.4964),
                (6.4470, -3.1069, -5.2207),
                (5.1393, -3.8409, -6.1385),
                (5.2166, -4.7142, -3.7828),
                (3.6966, -3.8914, -4.1361),
                (6.1042, -2.6245, -2.7671),
                (4.6823, -0.7125, -2.2719),
                (4.5102, 0.1523, -4.6210),
                (6.0517, -0.6121, -4.2723),
            )
        )
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