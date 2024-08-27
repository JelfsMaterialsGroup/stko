import numpy as np

import stko
from tests.molecular.subgroup.case_data import CaseData


def test_get_binder_distance(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.SubgroupAnalyser()

    results = analyser.calculate(case_data.molecule)
    for test in results:
        assert len(results[test]) == len(case_data.sub_group_data[test])
        assert np.allclose(
            results[test],
            case_data.sub_group_data[test],
            rtol=0,
            atol=1e-3,
        )
