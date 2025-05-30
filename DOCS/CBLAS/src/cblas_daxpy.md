# `cblas_daxpy.c` - C Interface to BLAS DAXPY routine

## Overview

This file provides a C interface (`cblas_daxpy`) to the Fortran BLAS Level 1 routine `daxpy`. The `daxpy` routine computes a scalar-vector product and adds the result to another vector: `Y = alpha*X + Y`. This C interface handles potential differences in integer types between C and Fortran and calls the underlying Fortran implementation.

## Key Components

### 1. `void API_SUFFIX(cblas_daxpy)(const CBLAS_INT N, const double alpha, const double *X, const CBLAS_INT incX, double *Y, const CBLAS_INT incY)`
   - **Description:** C wrapper for the `daxpy` BLAS routine. It performs the operation `Y = alpha*X + Y`.
   - **Parameters/Arguments:**
     - `N`: (`CBLAS_INT`) The number of elements in vectors `X` and `Y`.
     - `alpha`: (`double`) The scalar multiplier `alpha`.
     - `X`: (`const double *`) Pointer to the input vector `X`.
     - `incX`: (`CBLAS_INT`) The increment for the elements of vector `X`.
     - `Y`: (`double *`) Pointer to the input vector `Y`, which is overwritten with the result.
     - `incY`: (`CBLAS_INT`) The increment for the elements of vector `Y`.
   - **Returns:** `void`. The result is stored in the array pointed to by `Y`.

## Important Variables/Constants

- `F77_INT`: A macro used for type casting `CBLAS_INT` to the Fortran integer type if they differ (controlled by preprocessor directive `F77_INT`). If not defined, C integers are assumed to be compatible.
- `F77_N`, `F77_incX`, `F77_incY`: Local variables used to hold the (potentially type-casted) integer arguments passed to the Fortran `daxpy` routine.

## Usage Examples

**Example 1: Using `cblas_daxpy` in C**
```c
#include "cblas.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    CBLAS_INT n = 5;
    double alpha = 2.0;
    double *x = NULL;
    double *y = NULL;
    CBLAS_INT incx = 1;
    CBLAS_INT incy = 1;
    int i;

    x = (double*)malloc(n * sizeof(double));
    y = (double*)malloc(n * sizeof(double));

    if (!x || !y) {
        printf("Memory allocation failed\n");
        free(x);
        free(y);
        return 1;
    }

    // Initialize vectors
    // x = [1.0, 2.0, 3.0, 4.0, 5.0]
    // y = [10.0, 20.0, 30.0, 40.0, 50.0]
    for (i = 0; i < n; i++) {
        x[i] = (double)(i + 1);
        y[i] = (double)((i + 1) * 10);
    }

    printf("Vector X before cblas_daxpy:\n");
    for (i = 0; i < n; i++) printf("%.1f ", x[i]);
    printf("\n");

    printf("Vector Y before cblas_daxpy:\n");
    for (i = 0; i < n; i++) printf("%.1f ", y[i]);
    printf("\n");

    // Y = alpha*X + Y
    // Y = 2.0 * [1,2,3,4,5] + [10,20,30,40,50]
    // Y = [2,4,6,8,10] + [10,20,30,40,50]
    // Y = [12,24,36,48,60]
    API_SUFFIX(cblas_daxpy)(n, alpha, x, incx, y, incy);

    printf("Vector Y after cblas_daxpy:\n");
    for (i = 0; i < n; i++) printf("%.1f ", y[i]);
    printf("\n");

    free(x);
    free(y);

    return 0;
}

/*
To compile (assuming CBLAS library is linked, e.g., -lcblas):
gcc example_cblas_daxpy.c -o example_cblas_daxpy -lcblas
*/
```

## Dependencies and Interactions

- **Internal Dependencies:**
  - None within this specific C file beyond standard C features.
- **External Libraries/APIs:**
  - `cblas.h`: Contains the C declarations for CBLAS routines and types like `CBLAS_INT`.
  - `cblas_f77.h`: Contains prototypes for the Fortran BLAS routines callable from C (e.g., `F77_daxpy`).
  - The actual Fortran `daxpy` routine (typically from a BLAS library).
- **Called By:**
  - C applications that need to perform the `daxpy` operation using the CBLAS interface.
- **Calls To:**
  - `F77_daxpy`: The Fortran implementation of the `daxpy` routine.

---
*This documentation describes the C interface layer. For details on the underlying computation, refer to the documentation for the `daxpy` Fortran BLAS routine. Please refer to `CBLAS/src/cblas_daxpy.c` for the most up-to-date C interface code.*
```
