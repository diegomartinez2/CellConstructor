from __future__ import print_function
from __future__ import division

import numpy as np

import cellconstructor as CC
import cellconstructor.Settings
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.Units
import cellconstructor.symmetries

import pytest
import sys, os

try:
    import mpi4py.MPI
    _HAS_MPI = True
except ImportError:
    _HAS_MPI = False


def _get_comm():
    if _HAS_MPI:
        return mpi4py.MPI.COMM_WORLD
    return None


def _rank():
    return CC.Settings.get_rank()


def _nproc():
    return CC.Settings.GetNProc()


def _allclose_across_ranks(arr):
    """
    Return True if arr is identical on every MPI rank.
    MUST be called by ALL ranks simultaneously.
    """
    n = _nproc()
    if n == 1:
        return True
    comm = _get_comm()
    hashes = comm.gather(arr.ravel().tobytes()[:32], root=0)
    ok = True
    if _rank() == 0:
        ref = hashes[0]
        for h in hashes[1:]:
            if ref != h:
                ok = False
    ok = comm.bcast(ok, root=0)
    return ok


class TestTensorBroadcast:

    def test_broadcast_ndarray(self):
        """
        Master creates a numpy array, broadcasts it via
        Settings.broadcast(enforce_double=True). All ranks verify
        they received identical data.
        """
        shape = (27, 6, 6)
        if _rank() == 0:
            rng = np.random.RandomState(123)
            data = rng.randn(*shape).astype(np.float64)
        else:
            data = np.zeros(1)

        result = CC.Settings.broadcast(data, enforce_double=True)

        assert result.shape == shape
        assert result.dtype == np.float64
        assert np.max(np.abs(result)) > 0, "Broadcast produced all zeros"

        if _rank() == 0:
            assert np.allclose(result, data)

        assert _allclose_across_ranks(result), \
            "Broadcast data differs across MPI ranks!"

    def test_broadcast_ndarray_medium(self):
        """
        Broadcast a medium-sized array that exercises the broadcast
        path with a realistic data volume.
        """
        shape = (500, 100, 1)  # 50K float64 elements
        if _rank() == 0:
            rng = np.random.RandomState(456)
            data = rng.randn(*shape).astype(np.float64)
        else:
            data = np.zeros(1)

        result = CC.Settings.broadcast(data, enforce_double=True)

        assert result.shape == shape
        assert result.dtype == np.float64
        assert np.max(np.abs(result)) > 0

        if _rank() == 0:
            assert np.allclose(result, data)

        assert _allclose_across_ranks(result), \
            "Broadcast data differs across MPI ranks!"

    def _build_valid_dynmat(self, q_list):
        """
        Build dynamical matrices that satisfy D(-q)=D(q) so the
        real-space force constant is purely real.
        """
        nq = len(q_list)
        nat3 = 3 * 2
        rng = np.random.RandomState(42)

        dm = np.zeros((nq, nat3, nat3), dtype=np.complex128)
        paired = [False] * nq

        for i, qi in enumerate(q_list):
            if paired[i]:
                continue
            a = rng.randn(nat3, nat3)
            mat_real = 10.0 * 0.5 * (a + a.T)

            partner = None
            for j, qj in enumerate(q_list):
                if j == i:
                    continue
                if np.linalg.norm(np.array(qi) + np.array(qj)) < 1e-10:
                    partner = j
                    break

            if partner is not None and not paired[partner]:
                dm[i] = mat_real
                dm[partner] = mat_real
                paired[i] = paired[partner] = True
            else:
                dm[i] = mat_real
                paired[i] = True

        return dm

    def test_tensor2_setup_and_center(self):
        """
        Production-mimic test: builds a Tensor2 from a realistic Si
        supercell, calls SetupFromPhonons (broadcast #1) and Center
        (broadcast #2).  Verifies data is consistent across all ranks.
        """
        s = CC.Structure.Structure(2)
        s.atoms = ["Si", "Si"]
        s.coords[0, :] = [0.0,  0.0,  0.0]
        s.coords[1, :] = [0.25, 0.25, 0.25]
        s.unit_cell = np.eye(3) * 5.43
        s.has_unit_cell = True
        s.masses = {"Si": 28.0855}

        sc = (3, 3, 3)
        sc_struct = s.generate_supercell(sc)

        q_tot = CC.symmetries.GetQGrid(s.unit_cell, sc)
        nq = len(q_tot)

        phon = CC.Phonons.Phonons(s, nqirr=nq)
        phon.q_tot = [q.copy() for q in q_tot]
        phon.q_stars = [[q.copy() for q in q_tot]]

        dm = self._build_valid_dynmat(q_tot)
        for iq in range(nq):
            phon.dynmats[iq] = dm[iq]

        tensor = CC.ForceTensor.Tensor2(s, sc_struct, phon.GetSupercell())
        tensor.SetupFromPhonons(phon)

        before = tensor.tensor.copy()
        assert np.max(np.abs(before)) > 0

        assert _allclose_across_ranks(tensor.tensor), \
            "Tensor differs after SetupFromPhonons"

        tensor.Center()

        assert tensor.tensor.shape == before.shape
        assert np.max(np.abs(tensor.tensor)) > 0, \
            "Centered tensor is all zeros"
        assert tensor.x_r_vector2.shape[1] == tensor.n_R
        assert tensor.r_vector2.shape[1] == tensor.n_R

        assert _allclose_across_ranks(tensor.tensor), \
            "Tensor differs after Center"
        assert _allclose_across_ranks(tensor.x_r_vector2), \
            "x_r_vector2 differs after Center"
        assert _allclose_across_ranks(tensor.r_vector2), \
            "r_vector2 differs after Center"


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--pytest", action="store_true",
                        help="Run via pytest instead of standalone")
    parser.add_argument("-k", "--keyword", default=None,
                        help="pytest -k filter")
    args = parser.parse_args()

    if args.pytest or args.keyword:
        argv = [__file__, "-v"]
        if args.keyword:
            argv += ["-k", args.keyword]
        pytest.main(argv)
    else:
        runner = TestTensorBroadcast()
        r = _rank()
        n = _nproc()
        print("[rank {}/{}] === test_broadcast_ndarray ===".format(r, n))
        runner.test_broadcast_ndarray()
        print("[rank {}/{}] === test_broadcast_ndarray_medium ===".format(r, n))
        runner.test_broadcast_ndarray_medium()
        print("[rank {}/{}] === test_tensor2_setup_and_center ===".format(r, n))
        runner.test_tensor2_setup_and_center()
        print("[rank {}/{}] All tests passed.".format(r, n))
