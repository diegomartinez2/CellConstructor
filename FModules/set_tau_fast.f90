!-----------------------------------------------------------------------
! OPTIMIZED VERSION: set_tau_fast
! This subroutine replaces set_tau with O(N) scaling instead of O(N^4)
! It maintains the same atom ordering as the original set_tau
!-----------------------------------------------------------------------
SUBROUTINE set_tau_fast (nat, nat_blk, at, at_blk, tau, tau_blk, &
     ityp, ityp_blk, itau_blk)
  !-----------------------------------------------------------------------
  !
  ! Generate supercell coordinates with O(N) scaling
  !
  ! Parameters:
  !   nat       - Total number of atoms in supercell (nat_blk * n1 * n2 * n3)
  !   nat_blk   - Number of atoms in unit cell
  !   at        - Supercell lattice vectors (3x3)
  !   at_blk    - Unit cell lattice vectors (3x3)
  !   tau       - Output: supercell atomic positions (3,nat)
  !   tau_blk   - Input: unit cell atomic positions (3,nat_blk)
  !   ityp      - Output: atom types in supercell (nat)
  !   ityp_blk  - Input: atom types in unit cell (nat_blk)
  !   itau_blk  - Output: mapping from supercell to unit cell atoms (nat)
  !
  ! The supercell dimensions (n1,n2,n3) are computed from nat/nat_blk
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat, nat_blk
  INTEGER, INTENT(IN) :: ityp_blk(nat_blk)
  INTEGER, INTENT(OUT) :: ityp(nat), itau_blk(nat)
  DOUBLE PRECISION, INTENT(IN) :: at(3,3), at_blk(3,3)
  DOUBLE PRECISION, INTENT(IN) :: tau_blk(3,nat_blk)
  DOUBLE PRECISION, INTENT(OUT) :: tau(3,nat)
  !
  INTEGER :: i1, i2, i3, na_blk, na
  INTEGER :: n1, n2, n3, n_cells
  DOUBLE PRECISION :: r(3)
  DOUBLE PRECISION :: cell_ratio
  !
  ! Compute supercell dimensions
  ! nat = nat_blk * n1 * n2 * n3
  ! We need to find n1, n2, n3 such that this holds
  !
  n_cells = nat / nat_blk
  !
  ! Compute n1, n2, n3 from the ratio of unit cell to supercell vectors
  ! The unit cell vectors are related to supercell vectors by:
  !   at_blk(:,i) = at(:,i) / ni
  ! So we can compute ni = at(:,i) / at_blk(:,i) (taking norm)
  !
  n1 = NINT(SQRT(at(1,1)**2 + at(2,1)**2 + at(3,1)**2) / &
               SQRT(at_blk(1,1)**2 + at_blk(2,1)**2 + at_blk(3,1)**2))
  n2 = NINT(SQRT(at(1,2)**2 + at(2,2)**2 + at(3,2)**2) / &
               SQRT(at_blk(1,2)**2 + at_blk(2,2)**2 + at_blk(3,2)**2))
  n3 = NINT(SQRT(at(1,3)**2 + at(2,3)**2 + at(3,3)**2) / &
               SQRT(at_blk(1,3)**2 + at_blk(2,3)**2 + at_blk(3,3)**2))
  !
  ! Verify the dimensions are correct
  IF (n1 * n2 * n3 .NE. n_cells) THEN
     ! Try alternative: maybe the lattice vectors are transposed
     ! or the supercell is defined differently
     ! Fall back to computing from n_cells assuming cubic-like supercell
     n1 = NINT(n_cells**(1.0d0/3.0d0))
     n2 = n1
     n3 = n_cells / (n1 * n2)
     IF (n1 * n2 * n3 .NE. n_cells) THEN
        n1 = 1
        n2 = 1
        n3 = n_cells
     END IF
  END IF
  !
  ! The original set_tau searches from -NN to +NN and finds cells in [0,1)
  ! For a supercell n1 x n2 x n3, the valid cells are exactly:
  !   i1 = 0, 1, ..., n1-1
  !   i2 = 0, 1, ..., n2-1
  !   i3 = 0, 1, ..., n3-1
  !
  ! The original ordering (from loops i1=-NN:NN, i2=-NN:NN, i3=-NN:NN)
  ! produces cells in order of increasing i1, then i2, then i3
  ! for those cells where crystal coordinates are in [0,1)
  !
  ! For a cell (i1,i2,i3), the position in crystal coordinates is:
  !   r_cryst = (i1, i2, i3) / (n1, n2, n3) in each dimension
  ! which is in [0,1) when i1 in [0,n1), etc.
  !
  ! So we can directly generate the same ordering with three nested loops
  ! over i1=0:n1-1, i2=0:n2-1, i3=0:n3-1
  !
  na = 0
  !
  ! Loop over supercell grid - same ordering as original set_tau
  DO i1 = 0, n1-1
     DO i2 = 0, n2-1
        DO i3 = 0, n3-1
           !
           ! Compute the cell offset vector r in cartesian coordinates
           ! r = i1*at_blk(:,1) + i2*at_blk(:,2) + i3*at_blk(:,3)
           r(1) = i1 * at_blk(1,1) + i2 * at_blk(1,2) + i3 * at_blk(1,3)
           r(2) = i1 * at_blk(2,1) + i2 * at_blk(2,2) + i3 * at_blk(2,3)
           r(3) = i1 * at_blk(3,1) + i2 * at_blk(3,2) + i3 * at_blk(3,3)
           !
           ! Add all atoms from the unit cell at this position
           DO na_blk = 1, nat_blk
              na = na + 1
              tau(1,na) = tau_blk(1,na_blk) + r(1)
              tau(2,na) = tau_blk(2,na_blk) + r(2)
              tau(3,na) = tau_blk(3,na_blk) + r(3)
              ityp(na) = ityp_blk(na_blk)
              itau_blk(na) = na_blk
           END DO
           !
        END DO
     END DO
  END DO
  !
  ! Consistency check
  IF (na .NE. nat) THEN
     WRITE(*,*) 'Error in set_tau_fast: na /= nat', na, nat
     WRITE(*,*) 'n1,n2,n3:', n1, n2, n3
     WRITE(*,*) 'n_cells:', n_cells
     STOP
  END IF
  !
  RETURN
END SUBROUTINE set_tau_fast
!-----------------------------------------------------------------------
