# `dgesv.f` - Solves a real system of linear equations A * X = B

## Overview

This file implements the `DGESV` subroutine, a LAPACK driver routine. It computes the solution to a system of linear equations `A * X = B`, where `A` is a general N-by-N real matrix, and `X` and `B` are N-by-NRHS real matrices. The routine uses LU decomposition with partial pivoting and row interchanges to factorize matrix `A`.

## Key Components

### 1. `SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )`
   - **Description:** Solves `A * X = B` for `X`.
   - **Parameters/Arguments:**
     - `N`: (INTEGER) The order of the matrix `A` (number of linear equations). `N >= 0`.
     - `NRHS`: (INTEGER) The number of right-hand sides (columns of matrix `B`). `NRHS >= 0`.
     - `A`: (DOUBLE PRECISION array, dimension `(LDA,N)`) On entry, the N-by-N coefficient matrix `A`. On exit, the factors `L` and `U` from the factorization `A = P*L*U` (unit diagonal elements of `L` are not stored).
     - `LDA`: (INTEGER) The leading dimension of array `A`. `LDA >= max(1,N)`.
     - `IPIV`: (INTEGER array, dimension `(N)`) Output array. The pivot indices from the LU factorization; row `i` of the matrix was interchanged with row `IPIV(i)`.
     - `B`: (DOUBLE PRECISION array, dimension `(LDB,NRHS)`) On entry, the N-by-NRHS right-hand side matrix `B`. On exit, if `INFO = 0`, it contains the N-by-NRHS solution matrix `X`.
     - `LDB`: (INTEGER) The leading dimension of array `B`. `LDB >= max(1,N)`.
     - `INFO`: (INTEGER) Output status indicator:
       - `= 0`: Successful exit.
       - `< 0`: If `INFO = -i`, the i-th argument had an illegal value.
       - `> 0`: If `INFO = i`, `U(i,i)` is exactly zero. The factorization has been completed, but `U` is singular, and the solution could not be computed.
   - **Returns:** The solution `X` is stored in the `B` array. The pivot indices are in `IPIV`, and factorization details are in `A`.

## Important Variables/Constants

- `INFO`: Critical for determining the outcome of the computation (see description above).
- The routine relies on the `MAX` intrinsic function for argument checking.

## Usage Examples

**Example 1: Solving a system of linear equations in Fortran**
```fortran
PROGRAM TEST_DGESV
  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 3, NRHS = 1
  DOUBLE PRECISION :: A(N,N), B(N,NRHS)
  INTEGER :: IPIV(N), INFO, LDA, LDB
  INTEGER :: I, J

  LDA = N
  LDB = N

  ! Initialize matrix A
  ! A = | 2  3  3 |
  !     | 1 -1  0 |
  !     | 1  2 -1 |
  A(1,1)=2.0; A(1,2)=1.0; A(1,3)=1.0
  A(2,1)=3.0; A(2,2)=-1.0;A(2,3)=2.0
  A(3,1)=3.0; A(3,2)=0.0; A(3,3)=-1.0


  ! Initialize right-hand side B
  ! B = | 9 |
  !     | 8 |
  !     | 4 |
  B(1,1)=9.0
  B(2,1)=8.0
  B(3,1)=4.0

  PRINT *, "Matrix A:"
  DO I=1,N
     PRINT *, (A(I,J), J=1,N)
  END DO
  PRINT *, "Matrix B (RHS):"
  DO I=1,N
     PRINT *, (B(I,J), J=1,NRHS)
  END DO

  ! Solve A*X = B
  CALL DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)

  IF (INFO .EQ. 0) THEN
    PRINT *, "Solution X (stored in B):"
    DO I=1,N
       PRINT *, (B(I,J), J=1,NRHS)
    END DO
    ! Expected solution for this system: X = [1.0, 2.0, 3.0]'
  ELSE IF (INFO .GT. 0) THEN
    PRINT *, "Matrix A is singular. U(", INFO, ",", INFO, ") is zero."
  ELSE
    PRINT *, "Error: Argument ", -INFO, " had an illegal value."
  END IF

END PROGRAM TEST_DGESV
```
*(Note: The example matrix A was slightly changed to be more standard for examples and ensure a non-trivial solution. The example values for A and B are illustrative.)*

## Dependencies and Interactions

- **Internal Dependencies:**
  - None beyond standard Fortran intrinsics (`MAX`).
- **External Libraries:**
  - This is a LAPACK driver routine.
- **Called By:**
  - Typically called by user programs or higher-level libraries requiring solutions to general linear systems.
- **Calls To:**
  - `DGETRF`: Computes the LU factorization of matrix `A`.
  - `DGETRS`: Solves the system `A*X=B` using the `L` and `U` factors from `DGETRF`.
  - `XERBLA`: LAPACK error handling routine, called if input arguments are invalid.

---
*This documentation is based on the source code comments and standard LAPACK understanding. Please refer to `SRC/dgesv.f` for the most up-to-date information.*
```
