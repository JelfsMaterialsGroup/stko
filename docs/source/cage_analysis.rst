Cage analysis
=============

Here, we highlight useful analyses for ``Cage`` molecules, with some
internal pore. Many of these classes or methods have been applied in our
published works in slightly different forms, but generalised here. The
classes are not special, in that they can apply to any ``stk.Molecule``,
like all other methods in :mod:`stko`.

.. toctree::
  :maxdepth: 1

  ThreeSiteFG <_autosummary/stko.functional_groups.ThreeSiteFG>
  ThreeSiteFactory <_autosummary/stko.functional_groups.ThreeSiteFactory>
  CNCFactory <_autosummary/stko.functional_groups.CNCFactory>
  CNNFactory <_autosummary/stko.functional_groups.CNNFactory>
  NNNFactory <_autosummary/stko.functional_groups.NNNFactory>
  DecomposeMOC <_autosummary/stko.molecule_analysis.DecomposeMOC>
  DitopicThreeSiteAnalyser <_autosummary/stko.molecule_analysis.DitopicThreeSiteAnalyser>
  ConstructedAnalyser <_autosummary/stko.molecule_analysis.ConstructedAnalyser>
  GeometryAnalyser <_autosummary/stko.molecule_analysis.GeometryAnalyser>


There is an example script with some of the analyses
`here <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/cage_analysis_example.py>`_
that we run through below.

Define your building blocks, using standard ``stk`` behaviour.
And then you can build a Pd-based cage. Optimised with ``stk.MCHammer()``.
This example structure is not great, but it does the job here. Imagine
you would do a better optimisation with one of our
`Optimizer <optimizers.html>`_ classes!

.. code-block:: python

    import stk
    import stko

    pd = stk.BuildingBlock(
        smiles="[Pd+2]",
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)
        ),
        position_matrix=[[0.0, 0.0, 0.0]],
    )
    ditopic_bb = stk.BuildingBlock(
        smiles="C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]", bonders=(1,), deleters=()
            )
        ],
    )

    apdcage = stk.ConstructedMolecule(
        topology_graph=stk.cage.M6L12Cube(
            building_blocks=(pd, ditopic_bb),
            reaction_factory=stk.DativeReactionFactory(
                stk.GenericReactionFactory(
                    bond_orders={
                        frozenset(
                            {
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom,
                            }
                        ): 9,
                    },
                ),
            ),
            optimizer=stk.MCHammer(num_steps=1500),
        ),
    )


.. moldoc::

    import moldoc.molecule as molecule
    import stk

    pd = stk.BuildingBlock(
        smiles="[Pd+2]",
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)
        ),
        position_matrix=[[0.0, 0.0, 0.0]],
    )
    ditopic_bb = stk.BuildingBlock(
        smiles="C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]", bonders=(1,), deleters=()
            )
        ],
    )

    apdcage = stk.ConstructedMolecule(
        topology_graph=stk.cage.M6L12Cube(
            building_blocks=(pd, ditopic_bb),
            reaction_factory=stk.DativeReactionFactory(
                stk.GenericReactionFactory(
                    bond_orders={
                        # Use bond order of 1 here so that the
                        # rendering does not show a bond order
                        # of 9.
                        frozenset({
                            stk.GenericFunctionalGroup,
                            stk.SingleAtom,
                        }): 1,
                    },
                ),
            ),
            optimizer=stk.MCHammer(num_steps=1500),
        ),
    )

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            molecule.Atom(
                atomic_number=atom.get_atomic_number(),
                position=position,
            ) for atom, position in zip(
                apdcage.get_atoms(),
                apdcage.get_position_matrix(),
            )
        ),
        bonds=(
            molecule.Bond(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
                order=bond.get_order(),
            ) for bond in apdcage.get_bonds()
        ),
    )


We can decompose metal-organic cages (and any other really), by deleting
the metal atoms based on their atomic numbers. Thus producing a list of
distinct building blocks from the constructed molecule after optimisation.
Here, we use the ``stko.functional_groups.ThreeSiteFactory`` because the
functional group has three atoms, centered on a binder, which we can then
pass to the ``stko.molecule_analysis.DitopicThreeSiteAnalyser`` for
automatic analyses.

.. code-block:: python

    ligands = stko.molecule_analysis.DecomposeMOC().decompose(
        molecule=apdcage,
        metal_atom_nos=(46,),
    )
    print(f"there are {len(ligands)} ligands")

    tsa = stko.molecule_analysis.DitopicThreeSiteAnalyser()
    ligand_dict = defaultdict(list)
    for id_, lig in enumerate(ligands):
        # Write to file.
        lig.write(f"cage_output/apdcage_{id_}.mol")

        # Define as building block with
        as_building_block = stk.BuildingBlock.init_from_molecule(
            lig,
            stko.functional_groups.ThreeSiteFactory(smarts="[#6]~[#7X2]~[#6]"),
        )
        ligand_dict["binder_dist"].append(
            tsa.get_binder_distance(as_building_block)
        )
        ligand_dict["binder_binder"].append(
            tsa.get_binder_binder_angle(as_building_block)
        )
        ligand_dict["torsion"].append(
            abs(tsa.get_binder_adjacent_torsion(as_building_block))
        )
        binder_angles = tsa.get_binder_angles(as_building_block)
        ligand_dict["binder_angle"].append(binder_angles[0])
        ligand_dict["binder_angle"].append(binder_angles[1])
        ligand_dict["bite_angle"].append(
            sum(tsa.get_halfbite_angles(as_building_block))
        )

The printing the average of the collated values for all ligands looks
like:

.. code-block::

  there are 12 ligands
  avg. binder distance: 9.786611536300898
  avg. binder binder angle: 120.02568725689645
  avg. binder adjacent torsion: 0.9038734003287088
  avg. binder angles: 150.01204826771462
  avg. bite angle [caution]: 120.0240965354292

We can also analyse the cage itself, using some geometric measures. Note
that there are many external softwares for analysing cages too (like:
`pyWindow <https://github.com/JelfsMaterialsGroup/pywindow>`_). However,
we started with just simple measures here.

We can get the centroid and atom ids of distinct building blocks.
This simply extracts the parts of building blocks still present in the
molecule.

.. code-block:: python


    analyser = stko.molecule_analysis.ConstructedAnalyser()
    centroids = analyser.get_building_block_centroids(apdcage)
    # Out: {0: array([8.92788038, 0.41529785, 0.49671345]), ...}


We can get measures of pore size and cage geometry.

.. code-block:: python

    analyser = stko.molecule_analysis.GeometryAnalyser()
    print(
        f"approximate pore size: {analyser.get_min_centroid_distance(apdcage)}"
    )
    print(f"avg cage size: {analyser.get_avg_centroid_distance(apdcage)}")
    m_distances = list(
        analyser.get_metal_distances(
            apdcage,
            metal_atom_nos=(46,),
        ).values()
    )
    print(f"avg. metal distance: {np.mean(m_distances)}")
    m_angles = list(
        analyser.get_metal_centroid_metal_angle(
            apdcage,
            metal_atom_nos=(46,),
        ).values()
    )
    print(f"avg. metal-centroid-angles: {np.mean(m_angles)}")

    # And some geometrical measures.
    print(
        f"avg. N-Pd bond length:"
        f' {np.mean(analyser.calculate_bonds(apdcage)[("N", "Pd")])}'
    )
    print(
        f"N-Pd-N angles: "
        f'{analyser.calculate_angles(apdcage)[("N", "Pd", "N")]}'
    )

Giving:

.. code-block::

  approximate pore size: 6.905031217288509
  avg cage size: (9.796754899196744, 1.537005741468483)
  avg. metal distance: 13.760850522272712
  avg. metal-centroid-angles: 107.491677523429
  avg. N-Pd bond length: 2.0559406288465105
  N-Pd-N angles: [175.92946759351747, 89.4607655701432, 92.55953467090808, ...]
