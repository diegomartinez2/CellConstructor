import os

import numpy as np

import cellconstructor as CC
import cellconstructor.Phonons

# Two inequivalent atoms of the same species (Cu) in a tetragonal cell,
# 3x3x3 supercell. The reference dynamical matrices were produced by
# compute_phonons_finite_displacements driven by an exact harmonic force
# law built from this fixture's own FORCE_CONSTANTS (see
# generate_reference.py); that path has the known-correct atom assignment.
#
# Comparing the full matrices, not just the frequencies, is what makes this
# test sensitive to a wrong atom assignment: swapping the two atoms leaves
# every eigenvalue unchanged but shifts the matrix by ~1e-1, so a
# frequency-only check cannot see it.

HERE = os.path.dirname(os.path.abspath(__file__))


def _match_q(q_load, q_ref, unit_cell):
    """Index map so q_load[map[i]] equals q_ref[i] modulo a reciprocal vector."""
    bg = np.linalg.inv(unit_cell)
    mapping = []
    for q in q_ref:
        errs = [np.linalg.norm((d := (q - q2) @ bg) - np.round(d)) for q2 in q_load]
        j = int(np.argmin(errs))
        assert errs[j] < 1e-6, "q point %s not found in loaded grid" % (q,)
        mapping.append(j)
    return mapping


def test_load_phonopy_eigenvectors():
    dyn = CC.Phonons.Phonons()
    dyn.load_phonopy(os.path.join(HERE, "phonopy.yaml"))

    ref = np.load(os.path.join(HERE, "reference_dynmats.npy"))
    q_ref = np.load(os.path.join(HERE, "reference_qtot.npy"))
    mapping = _match_q(np.array(dyn.q_tot), q_ref, dyn.structure.unit_cell)

    max_diff = max(np.max(np.abs(ref[i] - dyn.dynmats[mapping[i]]))
                   for i in range(len(q_ref)))
    assert max_diff < 1e-3, \
        "load_phonopy dynamical matrix deviates from the reference by %.2e" % max_diff

    # The comparison must fail if the two atoms are swapped, otherwise it is
    # blind to the atom-assignment bug it is meant to catch.
    perm = np.arange(3 * dyn.structure.N_atoms)
    perm[0:3], perm[3:6] = [3, 4, 5], [0, 1, 2]
    swap_diff = max(np.max(np.abs(ref[i] - dyn.dynmats[mapping[i]][np.ix_(perm, perm)]))
                    for i in range(len(q_ref)))
    assert swap_diff > 1e-2, \
        "swapping the two atoms is not detected (%.2e); the test lacks teeth" % swap_diff


if __name__ == "__main__":
    test_load_phonopy_eigenvectors()
    print("ok")
