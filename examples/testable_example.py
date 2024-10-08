import os
from pathlib import Path

import aligner_example
import basic_example
import cage_analysis_example
import calculators_example
import mdanalysis_example
import molecule_splitter_example
import obabel_example
import openmm_example
import shape_example
import topology_extraction_example
import torsion_example
import zmatrix_example


def main() -> None:
    """Run the example."""
    init_dir = Path.cwd()
    os.chdir("examples/")

    try:
        aligner_example.main()
        basic_example.main()
        cage_analysis_example.main()
        calculators_example.main()
        mdanalysis_example.main()
        molecule_splitter_example.main()
        obabel_example.main()
        openmm_example.main()
        shape_example.main()
        topology_extraction_example.main()
        torsion_example.main()
        zmatrix_example.main()
    finally:
        os.chdir(init_dir)


if __name__ == "__main__":
    main()
