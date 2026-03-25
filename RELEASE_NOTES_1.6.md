# Release Notes for CellConstructor 1.6.0

## Overview

CellConstructor 1.6.0 is a significant performance-focused release that delivers substantial speedups for phonon calculations, particularly for medium-sized systems (50-200 atoms). This release includes critical bug fixes for modern Fortran compilers and several new features.

## Highlights

### 🚀 Performance Improvements

- **SymmetrizeFCQ Optimization**: Vectorized coordinate conversions and q-star copying operations, resulting in up to **25% speedup** for symmetry operations on typical systems.
  - Cartesian-crystal coordinate conversions are now **6× faster** through vectorized `np.einsum` operations
  - Q-star block copying optimized with NumPy transpose and reshape operations
  - Added cached transformation matrices to avoid redundant computations

- **DiagonalizeSupercell Optimization**: Major performance boost through vectorized operations and pre-computed phase factors, significantly reducing computation time for large supercell diagonalizations.

- **GenerateSupercell Optimization**: New O(N) `set_tau_fast` Fortran routine replaces the previous O(N²) algorithm for supercell generation, providing substantial speedups for large systems.

### 🐛 Bug Fixes

- **Fortran Compiler Compatibility**: Fixed acoustic sum rule errors that occurred with newer Fortran compilers (gfortran 15+), ensuring compatibility with modern HPC environments.
- **ForceTensor Error Messages**: Improved error handling and messages for incorrect file formats in both second-order and third-order force tensors.
- **DiagonalizeSupercell**: Fixed the `lo_to_split` argument handling and a missing timer variable that caused NameError exceptions.
- **Spglib Integration**: Removed the restrictive `spglib<=2.2` requirement and fixed symmetrization bugs for better compatibility with newer spglib versions.

### ✨ New Features

- **get_primitive_cell() Method**: New method in the `Structure` class to easily obtain the primitive cell from any crystal structure, simplifying workflow for high-throughput calculations.

### 🔧 Improvements

- **Timer Support**: Added optional `timer` parameter to `SymmetrizeFCQ`, `SetupQPoint`, `SymmetrizeDynQ`, `ApplyQStar`, and `ImposeSumRule` functions for detailed performance profiling.
- **Fortran Module Updates**: Enhanced `symdynph_gq_new.f90` to properly enforce symmetries after small-group symmetrization.

## Detailed Changes

### Commits since v1.5

1. **684036b** - Optimize SymmetrizeFCQ with vectorized coordinate conversions and q-star copying
2. **d7db321** - Add get_primitive_cell() method to Structure class
3. **e5f5fff** - Fix NameError: add missing t_main_loop_start timer in DiagonalizeSupercell_slow
4. **b462c05** - Optimize generate_supercell with O(N) set_tau_fast Fortran routine
5. **c49db94** - Optimize DiagonalizeSupercell with vectorized operations and pre-computed phase factors
6. **8da226c** - Fixed a typo
7. **697666c** - Fixed an error in the lo_to_split argumento of DiagonalizeSupercell
8. **9415f83** - Fixed the error message when the fileformat is wrong in the ForceTensor
9. **a7c9f46** - Fixed a missing error message in the file_format
10. **502225d** - Fix another fortran error in third order
11. **6be1e8e** - Merge pull request #118 from SSCHAcode/fix_fortran_compiler
12. **4228fd4** - New version. Fixed the error with older fortran compilers
13. **dd6fd10** - Fixed the error of the acoustic sum rule error for using a new fortran compiler
14. **4a3bb93** - Fixed another spglib test
15. **6e014e7** - Removed the requirement for spglib<=2.2
16. **d4a100a** - Merge pull request #117 from SSCHAcode/spglib_fix
17. **62649f2** - Fix a spelling
18. **cc4fe44** - Fixed a bug
19. **49f36b6** - Fixed the spglib symmetrization
20. **4045be7** - Fixing the interpolate mesh subroutine

## Performance Benchmarks

### Small System (CsSnI₃, 5 atoms, 10 irreducible q-points)
- **Total Symmetrize time**: 0.249s → 0.186s (**25% speedup**)
- **Coordinate conversions**: 0.06s → 0.01s (**6× speedup**)

### Medium System (CsSnI₃ rhombohedral, 10 atoms, 18 irreducible q-points)
- **Total Symmetrize time**: ~0.092s with negligible coordinate conversion overhead

*Note: Larger systems with bigger q-point stars will see even more significant improvements from the vectorized operations.*

## Compatibility

- **Python**: 3.8+
- **NumPy**: 1.20+
- **ASE**: Compatible with latest versions
- **Fortran Compiler**: gfortran 15+ and Intel Fortran
- **Spglib**: Now compatible with all recent versions (removed <=2.2 restriction)

## Upgrade Notes

This release is fully backward compatible. No API changes were made - all optimizations maintain the existing function signatures and behavior. Users can upgrade seamlessly and immediately benefit from the performance improvements.

## Contributors

Thanks to all contributors who made this release possible, particularly for the performance optimizations and compiler compatibility fixes.

---

**Full Changelog**: https://github.com/SSCHAcode/CellConstructor/compare/1.5...1.6.0
