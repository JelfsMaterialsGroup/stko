Using stk topology to facilitate GFN-FF optimisations
=====================================================

One of the key challenges with performing geometry optimisations from ``stk``
molecules is that they are often very far from equilibrium, meaning methods
like GFN2-xTB and GFN-FF, which are quite transferrable, will fail as they
infer the topology from the structure.

But ``stk`` molecules contain the bonding topology, so we can use that in
the optimisation, as we do here using the :class:`stko.GFNFFTopo` class.

.. testcode:: gfnfftopo-cage

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
        ),
    )

With no optimisation:

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

Then we can run the GFN-FF stk-topology informed optimisation:

.. testcode:: gfnfftopo-cage

    gfnff = stko.OptimizerSequence(
        stko.GFNFFTopo(
            version="harmonic2020",
            charge=12,
            solvent="dmso",
            square_planar=True,
        ),
        stko.GFNFFTopo(
            version="conformer2020",
            charge=12,
            solvent="dmso",
            square_planar=True,
        ),
    )
    apdcage = gfnff.optimize(apdcage)

And for the cost of about 4 seconds, you get this structure with GFN-FF
optimisation:

.. moldoc::

    import moldoc.molecule as molecule
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
        ),
    )
    gfnff = stko.OptimizerSequence(
        stko.GFNFFTopo(
            version="harmonic2020",
            charge=12,
            solvent="dmso",
            square_planar=True,
        ),
        stko.GFNFFTopo(
            version="conformer2020",
            charge=12,
            solvent="dmso",
            square_planar=True,
        ),
    )
    apdcage = gfnff.optimize(apdcage)

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
