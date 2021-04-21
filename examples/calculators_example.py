import stk
import stko


def main():
    bb1 = stk.BuildingBlock('NCCNCCN', [stk.PrimaryAminoFactory()])

    # Run optimisations.
    uff = stko.UFFEnergy()
    uff_results = uff.get_results(bb1)
    mmff = stko.MMFFEnergy()
    mmff_results = mmff.get_results(bb1)

    print(
        uff_results,
        uff_results.get_energy(),
        uff_results.get_unit_string(),
    )
    print(
        mmff_results,
        mmff_results.get_energy(),
        mmff_results.get_unit_string(),
    )


if __name__ == "__main__":
    main()
