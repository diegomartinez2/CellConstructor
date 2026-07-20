import cellconstructor as CC, cellconstructor.Phonons
import numpy as np
import sys, os


def test_consistency():
    """
    Verify that ExtractRandomStructures and ExtractRandomStructures_slow
    produce identical results.
    """
    test_dir = os.path.dirname(os.path.abspath(__file__))
    tdin_path = os.path.join(test_dir, "tdin")

    np.random.seed(42)
    dyn = CC.Phonons.Phonons(tdin_path, 40)
    np.random.seed(42)
    results_fast = dyn.ExtractRandomStructures(20)

    np.random.seed(42)
    dyn2 = CC.Phonons.Phonons(tdin_path, 40)
    np.random.seed(42)
    results_slow = dyn2.ExtractRandomStructures_slow(20)

    assert len(results_fast) == len(results_slow) == 20
    for i in range(20):
        max_diff = np.max(np.abs(results_fast[i].coords - results_slow[i].coords))
        assert max_diff < 1e-10, f"Structure {i}: coords differ (max diff = {max_diff})"


if __name__ == "__main__":
    test_consistency()
