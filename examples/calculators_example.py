import os
import sys

import stk
import stko


def main():
    if len(sys.argv) > 1:
        xtb_path = sys.argv[1]
    else:
        xtb_path = None

    bb1 = stk.BuildingBlock("NCCNCCN", [stk.PrimaryAminoFactory()])

    # Run calculations.
    calculations = []
    uff = stko.UFFEnergy()
    calculations.append(uff.calculate(bb1))
    uff_results = uff.get_results(bb1)
    mmff = stko.MMFFEnergy()
    calculations.append(mmff.calculate(bb1))
    mmff_results = mmff.get_results(bb1)

    print(
        uff_results,
        # Get energy from results.
        uff_results.get_energy(),
        uff_results.get_unit_string(),
        # Get energy directly through Calculator.
        uff.get_energy(bb1),
    )
    print(
        mmff_results,
        # Get energy from results.
        mmff_results.get_energy(),
        mmff_results.get_unit_string(),
        # Get energy directly through Calculator.
        mmff.get_energy(bb1),
    )

    if xtb_path is not None:
        print("doing XTB calculation.")
        xtb = stko.XTBEnergy(
            xtb_path=xtb_path,
            output_dir=os.path.join("output_directory", "example_xtb_out"),
            unlimited_memory=True,
            calculate_ip_and_ea=True,
        )

        xtb_results = xtb.get_results(bb1)
        print(xtb_results)
        # Extract properties from the energy calculator for a given
        # molecule.
        total_energy = xtb_results.get_total_energy()
        homo_lumo_gap = xtb_results.get_homo_lumo_gap()
        fermi_levels = xtb_results.get_fermi_level()
        homo_lumo_orbitals = xtb_results.get_homo_lumo_orbitals()
        full_dipole_moments = xtb_results.get_full_dipole_moments()
        ip = xtb_results.get_ionisation_potential()
        ea = xtb_results.get_electron_affinity()
        print(
            total_energy,
            homo_lumo_gap,
            homo_lumo_orbitals,
            fermi_levels,
            full_dipole_moments,
            ip,
            ea,
        )
        # From results, vs from calculator.
        print(xtb.get_energy(bb1), total_energy)
        try:
            xtb_results.get_total_free_energy()
        except AttributeError:
            print("Expected fail")

        # Try yielded option.
        calculations.append(xtb.calculate((bb1)))

    # Run through yield statements using the `.calculate` method
    # all in once.
    for i in calculations:
        print(f"{i}: xtb will not output anything as you need a new Results.")
        ey = stko.EnergyResults(i, "kcal mol-1")
        print(ey.get_energy(), ey.get_unit_string())

    # Zipping calculations together.
    mols = [
        stk.BuildingBlock("C1CCC1"),
        stk.BuildingBlock("C1CCCC1"),
        stk.BuildingBlock("C1CCCCC1"),
    ]
    for mol in mols:
        print(
            mol,
            stko.EnergyResults(mmff.calculate(mol), "kcal mol-1").get_energy(),
        )


if __name__ == "__main__":
    main()
