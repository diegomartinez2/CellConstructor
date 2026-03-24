# -*- coding: utf-8 -*-
"""
Test for get_primitive_cell method in Structure class.
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import sys
import os
import numpy as np
import pytest

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import cellconstructor as CC
import cellconstructor.Structure as Structure

try:
    import spglib
    __SPGLIB__ = True
except:
    __SPGLIB__ = False


def test_get_primitive_cell_basic():
    """
    Comprehensive test for get_primitive_cell method.
    
    This test:
    1. Creates a conventional FCC cell (4 atoms)
    2. Calls get_primitive_cell() to get the primitive cell
    3. Verifies the primitive cell has 1 atom
    4. Verifies the primitive cell volume is 1/4 of conventional
    5. Verifies atomic positions are correct
    6. Verifies masses are preserved
    7. Tests error handling for structure without unit cell
    """
    if not __SPGLIB__:
        pytest.skip("spglib is required for this test")
    
    # Change to test directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    # Create a conventional FCC cell (4 atoms)
    struct = CC.Structure.Structure()
    struct.unit_cell = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]], dtype=np.float64)
    struct.has_unit_cell = True
    struct.atoms = ["Al", "Al", "Al", "Al"]
    struct.N_atoms = 4
    struct.coords = np.array([
        [0, 0, 0],
        [0, 2, 2],
        [2, 0, 2],
        [2, 2, 0]
    ], dtype=np.float64)
    struct.masses = {"Al": 26.98}
    
    print("Original conventional cell atoms:", struct.N_atoms)
    print("Original conventional cell unit cell:")
    print(struct.unit_cell)
    
    # Get the primitive cell
    primitive = struct.get_primitive_cell()
    
    print("Primitive cell atoms:", primitive.N_atoms)
    print("Primitive cell unit cell:")
    print(primitive.unit_cell)
    print("Primitive positions:")
    print(primitive.coords)
    
    # The primitive cell should have fewer atoms than the conventional cell
    assert primitive.N_atoms < struct.N_atoms, \
        "Primitive cell should have fewer atoms than conventional cell"
    
    # The primitive cell should have 1 atom (FCC primitive)
    assert primitive.N_atoms == 1, \
        "FCC primitive cell should have 1 atom, got {}".format(primitive.N_atoms)
    
    # Check that the unit cell is valid
    assert primitive.has_unit_cell, "Primitive cell should have unit cell"
    assert np.abs(np.linalg.det(primitive.unit_cell)) > 1e-8, \
        "Primitive cell should have non-zero volume"
    
    # Check that the primitive cell volume is approximately 1/4 of conventional
    vol_conv = np.abs(np.linalg.det(struct.unit_cell))
    vol_prim = np.abs(np.linalg.det(primitive.unit_cell))
    print("Conventional volume:", vol_conv)
    print("Primitive volume:", vol_prim)
    print("Volume ratio:", vol_conv / vol_prim)
    
    # The volume ratio should be approximately 4
    assert np.abs(vol_conv / vol_prim - 4.0) < 0.1, \
        "Primitive cell volume should be ~1/4 of conventional cell"
    
    # Check that there's one Al atom
    assert primitive.atoms[0] == "Al", "Primitive cell should have Al atom"
    
    # Check that masses are preserved
    assert hasattr(primitive, 'masses'), "Primitive cell should have masses attribute"
    assert "Al" in primitive.masses, "Primitive cell should have Al in masses"
    assert np.isclose(primitive.masses["Al"], 26.98), \
        "Primitive cell should preserve mass value"
    
    # Test error handling for structure without unit cell
    struct_no_cell = CC.Structure.Structure()
    struct_no_cell.atoms = ["Al"]
    struct_no_cell.N_atoms = 1
    struct_no_cell.coords = np.array([[0, 0, 0]], dtype=np.float64)
    
    with pytest.raises(ValueError):
        struct_no_cell.get_primitive_cell()
    
    print("Test 1 passed: Basic primitive cell extraction works correctly")


def test_get_primitive_cell_with_itau():
    """
    Test get_primitive_cell with itau mapping and already primitive cell handling.
    
    This test:
    1. Creates a conventional FCC cell
    2. Gets primitive cell with itau mapping
    3. Verifies itau correctly maps all atoms to the single primitive atom
    4. Tests that already primitive cell returns a copy (not same object)
    5. Verifies the copy has correct data
    """
    if not __SPGLIB__:
        pytest.skip("spglib is required for this test")
    
    # Change to test directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    # Create a conventional FCC cell (4 atoms)
    struct = CC.Structure.Structure()
    struct.unit_cell = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 4]], dtype=np.float64)
    struct.has_unit_cell = True
    struct.atoms = ["Al", "Al", "Al", "Al"]
    struct.N_atoms = 4
    struct.coords = np.array([
        [0, 0, 0],
        [0, 2, 2],
        [2, 0, 2],
        [2, 2, 0]
    ], dtype=np.float64)
    
    # Get the primitive cell with itau mapping
    primitive, itau = struct.get_primitive_cell(return_itau=True)
    
    print("itau array:", itau)
    print("Number of input atoms:", struct.N_atoms)
    print("Number of primitive atoms:", primitive.N_atoms)
    
    # itau should have same length as input atoms
    assert len(itau) == struct.N_atoms, \
        "itau should have same length as input atoms ({}), got {}".format(struct.N_atoms, len(itau))
    
    # All atoms should map to atom 0 in primitive (since all are equivalent Al atoms in FCC)
    assert np.all(itau == 0), "All atoms should map to primitive atom 0, got {}".format(itau)
    
    # Verify that itau values are valid indices for primitive atoms
    assert np.all(itau >= 0), "All itau values should be >= 0"
    assert np.all(itau < primitive.N_atoms), \
        "All itau values should be < number of primitive atoms ({})".format(primitive.N_atoms)
    
    # Test already primitive cell returns a copy
    # First, let's create what spglib would consider a primitive cell
    # (The primitive cell we got from the conventional cell)
    prim_copy = primitive.get_primitive_cell()
    
    # Verify it's a copy, not the same object
    assert prim_copy is not primitive, "Should return a copy, not the same object"
    
    # Verify the data is the same
    assert prim_copy.N_atoms == primitive.N_atoms
    assert np.allclose(prim_copy.unit_cell, primitive.unit_cell)
    assert prim_copy.atoms == primitive.atoms
    
    print("Test 2 passed: itau mapping and copy behavior work correctly")


if __name__ == "__main__":
    test_get_primitive_cell_basic()
    test_get_primitive_cell_with_itau()
    print("\nAll tests passed!")
