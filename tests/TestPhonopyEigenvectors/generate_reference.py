"""Regenerate the phonopy fixture and the reference dynamical matrices.

Run once to (re)create phonopy.yaml, FORCE_CONSTANTS, reference_dynmats.npy
and reference_qtot.npy for test_phonopy_eigenvectors. Needs phonopy and ase;
the test itself needs neither.

Two inequivalent atoms of the same species (Cu) in a tetragonal cell make a
wrong atom assignment observable: the second atom sits off the mirror plane
so the two sites differ. A genuine phonopy FORCE_CONSTANTS/phonopy.yaml pair
is produced from EMT forces. An exact harmonic force law built from those
same force constants then drives compute_phonons_finite_displacements, whose
atom assignment is known correct; with a linear force law the finite
displacements recover the force constants exactly.
"""
import os

import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.calculators.calculator import Calculator, all_changes

import cellconstructor as CC
import cellconstructor.Phonons
from cellconstructor.Structure import Structure

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.file_IO import write_FORCE_CONSTANTS

HERE = os.path.dirname(os.path.abspath(__file__))
SC = (3, 3, 3)
EPS = 0.01
A, C = 3.6, 4.2
CELL = np.array([[A, 0, 0], [0, A, 0], [0, 0, C]], dtype=float)
SCALED = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.35]], dtype=float)
SYMBOLS = ["Cu", "Cu"]


def make_phonopy_fixture():
    unit = PhonopyAtoms(symbols=SYMBOLS, scaled_positions=SCALED, cell=CELL)
    phonon = Phonopy(unit, supercell_matrix=np.diag(SC))
    phonon.generate_displacements(distance=EPS)
    forces = []
    for sc in phonon.supercells_with_displacements:
        atoms = Atoms(symbols=sc.symbols, scaled_positions=sc.scaled_positions,
                      cell=sc.cell, pbc=True)
        atoms.calc = EMT()
        f = atoms.get_forces()
        forces.append(f - f.mean(axis=0))
    phonon.forces = forces
    phonon.produce_force_constants()
    phonon.save(os.path.join(HERE, "phonopy.yaml"),
                settings={"force_constants": False})
    # Write the compact format (one block row per primitive atom): it is far
    # smaller and forces load_phonopy through the reduced_to expansion, the
    # path where a wrong atom assignment can hide.
    p2s_map = phonon.primitive.p2s_map
    write_FORCE_CONSTANTS(phonon.force_constants[p2s_map],
                          os.path.join(HERE, "FORCE_CONSTANTS"),
                          p2s_map=p2s_map)
    return phonon


class HarmonicCalc(Calculator):
    """F = -Phi u, matching each atom to a reference site by minimum image so
    the result is independent of the supercell atom ordering."""
    implemented_properties = ["energy", "forces"]

    def __init__(self, ref_positions, ref_cell, fc_full, **kw):
        super().__init__(**kw)
        self.R0 = np.asarray(ref_positions, dtype=float)
        self.H = np.asarray(ref_cell, dtype=float)
        self.Hinv = np.linalg.inv(self.H)
        self.PHI = np.asarray(fc_full, dtype=float)
        self.n = len(self.R0)

    def calculate(self, atoms=None, properties=("energy",),
                  system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        idx = np.empty(self.n, dtype=int)
        U = np.zeros((self.n, 3))
        for k, p in enumerate(atoms.get_positions()):
            dmin = ((p[None, :] - self.R0) @ self.Hinv.T)
            dmin = (dmin - np.round(dmin)) @ self.H
            r = int(np.argmin(np.linalg.norm(dmin, axis=1)))
            idx[k], U[r] = r, dmin[r]
        forces_ref = -(self.PHI @ U.ravel()).reshape(self.n, 3)
        forces = np.array([forces_ref[r] for r in idx])
        self.results = {"energy": 0.5 * float(U.ravel() @ self.PHI @ U.ravel()),
                        "forces": forces}


def full_fc(phonon):
    fc = phonon.force_constants
    n = len(phonon.supercell)
    return fc.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)


def main():
    phonon = make_phonopy_fixture()
    calc = HarmonicCalc(phonon.supercell.positions, phonon.supercell.cell,
                        full_fc(phonon))
    struct = Structure()
    struct.generate_from_ase_atoms(
        Atoms(symbols=SYMBOLS, scaled_positions=SCALED, cell=CELL, pbc=True))
    dyn = CC.Phonons.compute_phonons_finite_displacements(
        struct, calc, epsilon=EPS, supercell=SC, use_symmetries=False)
    np.save(os.path.join(HERE, "reference_dynmats.npy"), np.array(dyn.dynmats))
    np.save(os.path.join(HERE, "reference_qtot.npy"), np.array(dyn.q_tot))
    print("regenerated fixture and reference for", len(dyn.q_tot), "q points")


if __name__ == "__main__":
    main()
