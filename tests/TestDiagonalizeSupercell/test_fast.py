"""
Test to verify DiagonalizeSupercell (new fast version) produces identical results to DiagonalizeSupercell_slow (original)
"""
from __future__ import print_function
import sys
import os
import numpy as np
import pytest

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.Timer



def test_diagonalize_supercell_fast():
    """
    Test that DiagonalizeSupercell (new fast version) matches DiagonalizeSupercell_slow (original) exactly.
    
    This test:
    1. Loads a dynamical matrix
    2. Runs both DiagonalizeSupercell (new fast) and DiagonalizeSupercell_slow (original)
    3. Compares frequencies (should match within numerical precision)
    4. Compares polarization vectors (allowing for sign flips which are physically equivalent)
    """
    # Go to the test directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    print("\n" + "="*60)
    print("Testing DiagonalizeSupercell (new) vs DiagonalizeSupercell_slow (original)")
    print("="*60)
    
    # Load the test dynamical matrix (CsSnI3 with 10 irreducible q-points)
    print("\nLoading dynamical matrix...")
    dyn = CC.Phonons.Phonons("CsSnI3_dyn_", 10)
    print(f"Loaded: {dyn.structure.N_atoms} atoms, {len(dyn.q_tot)} q-points")
    
    # Run the original (slow) version
    print("\nRunning DiagonalizeSupercell_slow (original)...")
    timer = CC.Timer.Timer(active=True)
    w_slow, p_slow = timer.execute_timed_function(dyn.DiagonalizeSupercell_slow)
    t_slow = timer.timed_subroutines['DiagonalizeSupercell_slow']['time']
    print(f"Original (slow) version took: {t_slow:.4f} seconds")
    
    # Reload to ensure clean state (avoid any modifications from first run)
    dyn2 = CC.Phonons.Phonons("CsSnI3_dyn_", 10)
    
    # Run the new (fast) version
    print("\nRunning DiagonalizeSupercell (new fast)...")
    timer2 = CC.Timer.Timer(active=True)
    w_fast, p_fast = timer2.execute_timed_function(dyn2.DiagonalizeSupercell)
    t_fast = timer2.timed_subroutines['DiagonalizeSupercell']['time']
    print(f"New (fast) version took: {t_fast:.4f} seconds")
    
    # Report speedup
    speedup = t_slow / t_fast if t_fast > 0 else float('inf')
    print(f"\nSpeedup: {speedup:.2f}x")
    
    # Compare frequencies
    print("\nComparing frequencies...")
    freq_diff = np.max(np.abs(w_slow - w_fast))
    print(f"Max frequency difference: {freq_diff:.2e}")
    
    # The frequencies should match very closely
    assert freq_diff < 1e-8, f"Frequency mismatch too large: {freq_diff}"
    
    # Compare polarization vectors
    print("\nComparing polarization vectors...")
    max_pol_diff = 0.0
    modes_with_diff = 0
    
    for i in range(p_slow.shape[1]):
        p_s = p_slow[:, i]
        p_f = p_fast[:, i]
        
        # Check direct match or sign-flipped match (both are physically valid)
        diff_direct = np.max(np.abs(p_s - p_f))
        diff_flipped = np.max(np.abs(p_s + p_f))
        min_diff = min(diff_direct, diff_flipped)
        
        if min_diff > max_pol_diff:
            max_pol_diff = min_diff
            
        if min_diff > 1e-8:
            modes_with_diff += 1
            if modes_with_diff <= 5:  # Print first 5 problematic modes
                print(f"  Mode {i}: direct_diff={diff_direct:.2e}, flipped_diff={diff_flipped:.2e}")
    
    print(f"Max polarization difference: {max_pol_diff:.2e}")
    print(f"Modes with significant difference: {modes_with_diff}")
    
    # The polarization vectors should match (allowing for sign flips)
    assert max_pol_diff < 1e-8, f"Polarization vector mismatch too large: {max_pol_diff}"
    
    print("\n" + "="*60)
    print("TEST PASSED: Both functions produce identical results!")
    print(f"Speedup achieved: {speedup:.2f}x")
    print("="*60 + "\n")


if __name__ == "__main__":
    test_diagonalize_supercell_fast()
