import os

import cellconstructor as CC
import cellconstructor.Phonons

import numpy as np

# Frequencies (cm^-1) of the bcc Nb fixture on the commensurate q grid,
# computed independently with phonopy 2.x from the same
# phonopy.yaml / FORCE_CONSTANTS pair.
REFERENCE_FREQS_CM = np.sort(np.concatenate([
    [0.0] * 3,
    [91.427] * 12,
    [138.125] * 6,
    [168.055] * 6,
    [172.073] * 6,
    [176.968] * 6,
    [195.129] * 6,
    [217.426] * 3,
]))

# Lattice parameter of the conventional cell, given as 6.19929787 au
# in the yaml file.
A_CONV_ANGSTROM = 6.199297870 * CC.Units.BOHR_TO_ANGSTROM


def check_dyn(dyn):
    assert np.allclose(dyn.structure.unit_cell, np.eye(3) * A_CONV_ANGSTROM,
                       atol=1e-6)
    assert np.allclose(dyn.structure.coords[1, :],
                       np.ones(3) * A_CONV_ANGSTROM / 2, atol=1e-6)

    w, pols = dyn.DiagonalizeSupercell()
    freqs_cm = np.sort(w * CC.Units.RY_TO_CM)
    assert np.allclose(freqs_cm, REFERENCE_FREQS_CM, atol=0.01)


def test_phonopy_input():
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    dyn = CC.Phonons.Phonons()
    dyn.load_phonopy("phonopy.yaml")

    check_dyn(dyn)

    dyn.save_qe("prova")


def test_phonopy_input_ev_angstrom(tmp_path):
    # Convert the qe-unit fixture (au, Ry/au^2) to the phonopy default
    # units (Angstrom, eV/Angstrom^2): the loaded dynamical matrix must
    # be identical.
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    in_lattice = False
    yaml_lines = []
    with open("phonopy.yaml", "r") as fp:
        for line in fp:
            stripped = line.strip()
            if stripped.endswith(":"):
                in_lattice = stripped == "lattice:"
            elif in_lattice and stripped.startswith("-"):
                values = [float(x) for x in
                          stripped.split("[")[1].split("]")[0].split(",")]
                values = [v * CC.Units.BOHR_TO_ANGSTROM for v in values]
                line = "  - [ %22.15f, %22.15f, %22.15f ]\n" % tuple(values)
            elif stripped.startswith("length:"):
                line = '  length: "angstrom"\n'
            elif stripped.startswith("force_constants:"):
                line = '  force_constants: "eV/angstrom^2"\n'
            yaml_lines.append(line)

    fc_lines = []
    fc_factor = CC.Units.RY_TO_EV / CC.Units.BOHR_TO_ANGSTROM**2
    with open("FORCE_CONSTANTS", "r") as fp:
        for line in fp:
            data = line.split()
            if len(data) == 3:
                values = [float(x) * fc_factor for x in data]
                line = " %20.12f %20.12f %20.12f\n" % tuple(values)
            fc_lines.append(line)

    yaml_filename = str(tmp_path / "phonopy.yaml")
    with open(yaml_filename, "w") as fp:
        fp.writelines(yaml_lines)
    with open(str(tmp_path / "FORCE_CONSTANTS"), "w") as fp:
        fp.writelines(fc_lines)

    dyn = CC.Phonons.Phonons()
    dyn.load_phonopy(yaml_filename)

    check_dyn(dyn)


if __name__ == "__main__":
    test_phonopy_input()
