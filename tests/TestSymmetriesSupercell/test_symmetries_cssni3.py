import cellconstructor as CC, cellconstructor.Phonons
import numpy as np
import sys, os


def test_symmetries_cssni3():
    os.chdir(os.path.dirname(__file__))

    # Load the dynamic matrix
    dyn = CC.Phonons.Phonons("CsSnI3_rhombo_", 18)

    w, p = dyn.DiagonalizeSupercell()
    dyn.Symmetrize()
    w2, p2 = dyn.DiagonalizeSupercell()

    for i in range(len(w)):
        assert np.allclose(w[i], w2[i], rtol=1e-5), f"Frequencies differ for mode {i}: {w[i]} vs {w2[i]}"

if __name__ == "__main__":
    test_symmetries_cssni3()
