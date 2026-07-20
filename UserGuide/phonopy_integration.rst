Phonopy Integration
********************

CellConstructor can import dynamical matrices from phonopy output files and
convert between phonopy and CellConstructor data structures.

.. contents::
   :local:
   :depth: 2


Loading from phonopy files
==========================

The method :meth:`Phonons.load_phonopy` loads the dynamical matrix directly
from the files produced by a phonopy force-constant calculation. This is the
recommended way when the phonopy output is already on disk (no need for the
phonopy Python library at runtime — only ``numpy`` is required).

Required files
--------------

Two files are needed:

- ``phonopy.yaml``
    Contains the unit cell, supercell, atom positions, masses, the
    supercell matrix, and optionally the physical units. The importer
    reads the ``unit_cell``, ``supercell``, ``supercell_matrix``, and
    ``physical_unit`` blocks.

- ``FORCE_CONSTANTS``
    The real-space force constant matrix in either the **full** format
    (one row per supercell atom) or the **compact** format (rows only for
    primitive-cell atoms, expanded via the ``reduced_to`` field of the
    YAML file).

Both files are written by ``phonopy`` when you save the results of a
force-constant calculation.

Basic usage
-----------

.. code:: python

   import cellconstructor as CC
   import cellconstructor.Phonons

   dyn = CC.Phonons.Phonons()
   dyn.load_phonopy("path/to/phonopy.yaml")

   # Now the dynamical matrix is loaded; use it normally:
   w, pols = dyn.DiagonalizeSupercell()

   # Save in Quantum ESPRESSO format:
   dyn.save_qe("dyn")

If the ``FORCE_CONSTANTS`` file is in a different location, pass the path
explicitly:

.. code:: python

   dyn.load_phonopy("path/to/phonopy.yaml",
                    fc_filename="other/path/FORCE_CONSTANTS")

Unit handling
-------------

The importer reads the ``physical_unit`` block from ``phonopy.yaml`` and
converts lengths to Angstrom and force constants to Ry/bohr\ :sup:`2`
(CellConstructor's internal conventions).

Supported length units:
  ``angstrom``, ``au``, ``bohr``

Supported force-constant units:
  ``eV/angstrom^2``, ``eV/angstrom.au``, ``eV/au^2``, ``Ry/au^2``,
  ``mRy/au^2``, ``hartree/au^2``

If the ``physical_unit`` block is missing, the phonopy defaults are assumed:

- Length: ``angstrom``
- Force constants: ``eV/angstrom^2``

For example, a VASP calculation produces phonopy files in Angstrom and
eV/Angstrom\ :sup:`2`, while a Quantum ESPRESSO calculation usually produces
files in Bohr (``au``) and Ry/au\ :sup:`2`. Both are handled automatically.

Full vs compact FORCE_CONSTANTS
-------------------------------

Phonopy can write ``FORCE_CONSTANTS`` in two formats:

- **Full format**: Each supercell atom has its own row in the file.
  The first header line reads ``N_sc  N_sc`` (where *N_sc* is the number of
  supercell atoms).

- **Compact format**: Only the primitive-cell atoms have rows.
  The first header line reads ``N_prim  N_sc`` (where *N_prim* is the number
  of primitive-cell atoms). The remaining rows are recovered using the
  translational symmetry

  .. math::

     \Phi(p + \mathbf{T}, \, j) = \Phi(p, \, j - \mathbf{T})

  where *p* is a primitive atom and *\mathbf{T}* the lattice translation
  that connects a supercell atom to its primitive representative (given by
  the ``reduced_to`` field in ``phonopy.yaml``).

Both formats are supported transparently — the importer detects the format
from the file header.

Limitations
-----------

- Only **diagonal** supercell matrices are supported (the common case for
  phonon calculations). Off-diagonal supercell matrices raise
  ``NotImplementedError``.
- The ``primitive_cell`` block in the YAML is not used; the
  ``unit_cell`` block must be present.
- The ``phonopy.yaml`` must contain the ``reduced_to`` fields when the
  compact ``FORCE_CONSTANTS`` format is used.


Loading from a phonopy Python object
=====================================

If you already have a ``phonopy.Phonopy`` object in memory, use
:func:`CC.Methods.sscha_phonons_from_phonopy`. This requires the phonopy
Python library at runtime (``import phonopy``).

.. code:: python

   import cellconstructor as CC
   from cellconstructor.Methods import sscha_phonons_from_phonopy

   # phonon is a phonopy.Phonopy object from your own code
   dyn = sscha_phonons_from_phonopy(phonon)

This returns a fully populated :class:`CC.Phonons.Phonons` object ready for
use.


Saving to phonopy format
=========================

.. note::

   The built-in method :meth:`Phonons.save_phonopy` is currently **broken**
   (raises ``NotImplementedError``). Use :func:`CC.Methods.tensor2_to_phonopy_fc2`
   and :func:`CC.Methods.tensor3_to_phonopy_fc3` instead.

Convert a second-order tensor to phonopy force constants:

.. code:: python

   from cellconstructor.Methods import tensor2_to_phonopy_fc2

   # SSCHA_tensor is a CC.ForceTensor.Tensor2 object
   # phonon is a phonopy.Phonopy object with the same supercell
   phonopy_fc2 = tensor2_to_phonopy_fc2(SSCHA_tensor, phonon)

The result is in eV/Angstrom\ :sup:`2` units, ready for use with phonopy or
phono3py.

Third-order force constants are exported similarly with
:func:`tensor3_to_phonopy_fc3 <CC.Methods.tensor3_to_phonopy_fc3>`, which
returns a tensor in eV/Angstrom\ :sup:`3`.

Both routines require the phonopy Python library.


Complete workflow example
=========================

Here is a full end-to-end example: load a phonopy calculation, compute
phonon frequencies, and export to Quantum ESPRESSO format.

.. code:: python

   import cellconstructor as CC
   import cellconstructor.Phonons

   # Load from phonopy
   dyn = CC.Phonons.Phonons()
   dyn.load_phonopy("phonopy.yaml")

   # Inspect the structure
   print("Number of atoms:", dyn.structure.N_atoms)
   print("Unit cell:", dyn.structure.unit_cell)

   # Diagonalize the supercell to get frequencies
   w, pols = dyn.DiagonalizeSupercell()

   # Convert from Rydberg to cm^-1
   freqs_cm = w * CC.Units.RY_TO_CM
   print("Frequencies (cm^-1):", freqs_cm)

   # Export to Quantum ESPRESSO format
   dyn.save_qe("dyn_qe")

   # The dynamical matrix can now be used for SSCHA,
   # Fourier interpolation, or any other CC feature.
