# `lapacke_dgesv.c` - C Interface to LAPACK DGESV routine

## Overview

This file provides a C interface (`LAPACKE_dgesv`) to the Fortran LAPACK driver routine `dgesv`. The `dgesv` routine computes the solution to a real system of linear equations `A * X = B`. This C interface handles matrix layout (row-major or column-major), performs optional NaN checks on input matrices, and then calls a working function (`LAPACKE_dgesv_work`) which further interfaces with the Fortran implementation.

## Key Components

### 1. `lapack_int API_SUFFIX(LAPACKE_dgesv)(int matrix_layout, lapack_int n, lapack_int nrhs, double* a, lapack_int lda, lapack_int* ipiv, double* b, lapack_int ldb)`
   - **Description:** C wrapper for the LAPACK `dgesv` routine. It solves `A * X = B` for `X`.
   - **Parameters/Arguments:**
     - `matrix_layout`: (int) Specifies whether matrices `a` and `b` are stored in row-major (`LAPACK_ROW_MAJOR`) or column-major (`LAPACK_COL_MAJOR`) order.
     - `n`: (`lapack_int`) The order of the matrix `A` (number of linear equations).
     - `nrhs`: (`lapack_int`) The number of right-hand sides (columns of matrix `B`).
     - `a`: (`double*`) Pointer to the input/output coefficient matrix `A`. On exit, contains `L` and `U` factors.
     - `lda`: (`lapack_int`) The leading dimension of array `a`.
     - `ipiv`: (`lapack_int*`) Pointer to the output array for pivot indices.
     - `b`: (`double*`) Pointer to the input/output right-hand side matrix `B`. On exit, if successful, contains the solution matrix `X`.
     - `ldb`: (`lapack_int`) The leading dimension of array `b`.
   - **Returns:** (`lapack_int`)
     - `0`: Successful exit.
     - `-1`: If `matrix_layout` is invalid.
     - `-4`: If matrix `a` contains NaN values (and NaN check is enabled).
     - `-7`: If matrix `b` contains NaN values (and NaN check is enabled).
     - Other negative values typically indicate an issue in the `_work` function or the underlying Fortran routine, often corresponding to the position of an invalid argument in the Fortran call.
     - Positive values: If `INFO = i` from the Fortran routine, `U(i,i)` is exactly zero, indicating singularity.

## Important Variables/Constants

- `LAPACK_COL_MAJOR`, `LAPACK_ROW_MAJOR`: Macros (typically defined in `lapacke.h` or `lapacke_utils.h`) used to specify matrix storage format.
- `LAPACK_DISABLE_NAN_CHECK`: Preprocessor macro that can be defined to disable NaN checking in LAPACKE functions.

## Usage Examples

**Example 1: Using `LAPACKE_dgesv` in C (column-major)**
```c
#include "lapacke.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    lapack_int n = 3, nrhs = 1;
    lapack_int lda = 3, ldb = 3;
    lapack_int info;
    lapack_int* ipiv = NULL;
    double* a = NULL;
    double* b = NULL;
    int i, j;

    // Allocate memory
    ipiv = (lapack_int*)malloc(n * sizeof(lapack_int));
    a = (double*)malloc(lda * n * sizeof(double));
    b = (double*)malloc(ldb * nrhs * sizeof(double));

    if (!ipiv || !a || !b) {
        printf("Memory allocation failed\n");
        free(ipiv); free(a); free(b);
        return 1;
    }

    // Initialize matrix A (column-major)
    // A = | 2  1  1 |
    //     | 3 -1  0 |
    //     | 3  2 -1 |
    // Transposed from typical row-major representation for column-major storage:
    // a = [2, 3, 3,  1, -1, 2,  1, 0, -1]
    a[0] = 2.0; a[1] = 3.0; a[2] = 3.0; // Column 1
    a[3] = 1.0; a[4] = -1.0;a[5] = 0.0; // Column 2
    a[6] = 1.0; a[7] = 2.0; a[8] = -1.0;// Column 3


    // Initialize right-hand side B (column-major)
    // B = | 9 |
    //     | 8 |
    //     | 4 |
    b[0] = 9.0; b[1] = 8.0; b[2] = 4.0;

    printf("Matrix A (column-major storage):\n");
    for(i=0; i<n; i++) { // Print row by row for easier reading
        for(j=0; j<n; j++) printf("%.1f ", a[j*lda+i]);
        printf("\n");
    }
    printf("Matrix B (RHS, column-major storage):\n");
    for(i=0; i<n; i++) {
         for(j=0; j<nrhs; j++) printf("%.1f ", b[j*ldb+i]);
         printf("\n");
    }


    // Solve A*X = B using column-major layout
    info = API_SUFFIX(LAPACKE_dgesv)(LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);

    if (info == 0) {
        printf("Solution X (stored in B, column-major):\n");
        for(i=0; i<n; i++) {
            for(j=0; j<nrhs; j++) printf("%.1f ", b[j*ldb+i]);
            printf("\n");
        }
        // For the A and B above, X should be approximately [1.0, 2.0, 3.0]'
    } else if (info > 0) {
        printf("Matrix A is singular. U(%d,%d) is zero.\n", (int)info, (int)info);
    } else {
        printf("Error: DGESV failed with info = %d\n", (int)info);
    }

    free(ipiv); free(a); free(b);
    return 0;
}

/*
To compile (assuming LAPACKE and BLAS libraries are linked, e.g., -llapacke -lcblas):
gcc example_lapacke_dgesv.c -o example_lapacke_dgesv -llapacke -lcblas -lm
*/
```
*(Note: The example matrix A was slightly changed for clarity in column-major representation.)*


## Dependencies and Interactions

- **Internal Dependencies:**
  - `lapacke_utils.h`: Provides utility functions and macros for LAPACKE, including `LAPACKE_xerbla` (error handler) and potentially `LAPACKE_dge_nancheck`.
- **External Libraries/APIs:**
  - The underlying Fortran `dgesv` LAPACK routine.
  - BLAS routines (implicitly used by `dgesv`).
- **Called By:**
  - C applications that need to solve systems of linear equations using the LAPACKE interface.
- **Calls To:**
  - `API_SUFFIX(LAPACKE_xerbla)`: For reporting errors related to invalid `matrix_layout`.
  - `API_SUFFIX(LAPACKE_dge_nancheck)`: (Optional, if NaN checking is enabled) To check for NaN values in input matrices `a` and `b`.
  - `API_SUFFIX(LAPACKE_dgesv_work)`: The workhorse function that handles data transformations (if needed for row-major) and calls the actual Fortran `dgesv` routine.

---
*This documentation describes the C interface layer. For details on the underlying computation, refer to the documentation for the `dgesv` Fortran LAPACK routine. Please refer to `LAPACKE/src/lapacke_dgesv.c` for the most up-to-date C interface code.*
```
