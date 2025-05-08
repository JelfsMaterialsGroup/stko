import contextlib
import logging
import os
from pathlib import Path

import aligner_example
import basic_example
import cage_analysis_example
import calculators_example
import intermediates_example
import mdanalysis_example
import molecule_splitter_example
import obabel_example
import openmm_example
import optwrite_example
import shape_example
import topology_extraction_example
import torsion_example
import zmatrix_example

from stko import WrapperNotInstalledError

logger = logging.getLogger(__name__)


def main() -> None:
    """Run the example."""
    init_dir = Path.cwd()
    os.chdir("examples/")

    try:
        intermediates_example.main()
        aligner_example.main()
        basic_example.main()
        cage_analysis_example.main()
        calculators_example.main()
        with contextlib.suppress(WrapperNotInstalledError):
            mdanalysis_example.main()
        molecule_splitter_example.main()
        with contextlib.suppress(WrapperNotInstalledError):
            obabel_example.main()
        openmm_example.main()
        shape_example.main()
        topology_extraction_example.main()
        torsion_example.main()
        zmatrix_example.main()
        optwrite_example.main()
        logger.info("all examples ran, at least!")

    finally:
        os.chdir(init_dir)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
