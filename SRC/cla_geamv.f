*> \brief \b CLA_GEAMV computes a matrix-vector product using a general matrix to calculate error bounds.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CLA_GEAMV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_geamv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_geamv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_geamv.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLA_GEAMV( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA,
*                             Y, INCY )
*
*       .. Scalar Arguments ..
*       REAL               ALPHA, BETA
*       INTEGER            INCX, INCY, LDA, M, N
*       INTEGER            TRANS
*       ..
*       .. Array Arguments ..
*       COMPLEX            A( LDA, * ), X( * )
*       REAL               Y( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLA_GEAMV  performs one of the matrix-vector operations
*>
*>         y := alpha*abs(A)*abs(x) + beta*abs(y),
*>    or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> m by n matrix.
*>
*> This function is primarily used in calculating error bounds.
*> To protect against underflow during evaluation, components in
*> the resulting vector are perturbed away from zero by (N+1)
*> times the underflow threshold.  To prevent unnecessarily large
*> errors for block-structure embedded in general matrices,
*> "symbolically" zero components are not perturbed.  A zero
*> entry is considered "symbolic" if all multiplications involved
*> in computing that entry have at least one zero multiplicand.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is INTEGER
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y)
*>             BLAS_TRANS         y := alpha*abs(A**T)*abs(x) + beta*abs(y)
*>             BLAS_CONJ_TRANS    y := alpha*abs(A**T)*abs(x) + beta*abs(y)
*>
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of the matrix A.
*>           M must be at least zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is REAL
*>           On entry, ALPHA specifies the scalar alpha.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,n)
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, m ).
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX array, dimension
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*>           Before entry, the incremented array X must contain the
*>           vector x.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is REAL
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*>           Unchanged on exit.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is REAL array, dimension
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*>           Before entry with BETA non-zero, the incremented array Y
*>           must contain the vector y. On exit, Y is overwritten by the
*>           updated vector y.
*>           If either m or n is zero, then Y not referenced and the function
*>           performs a quick return.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*>           Unchanged on exit.
*>
*>  Level 2 Blas routine.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup la_geamv
*
*  =====================================================================
      SUBROUTINE CLA_GEAMV( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                      BETA,
     $                      Y, INCY )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      INTEGER            TRANS
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * )
      REAL               Y( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            SYMB_ZERO
      REAL               TEMP, SAFE1
      INTEGER            I, INFO, IY, J, JX, KX, KY, LENX, LENY
      COMPLEX            CDUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, SLAMCH
      REAL               SLAMCH
*     ..
*     .. External Functions ..
      EXTERNAL           ILATRANS
      INTEGER            ILATRANS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, ABS, REAL, AIMAG, SIGN
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.( ( TRANS.EQ.ILATRANS( 'N' ) )
     $           .OR. ( TRANS.EQ.ILATRANS( 'T' ) )
     $           .OR. ( TRANS.EQ.ILATRANS( 'C' ) ) ) ) THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CLA_GEAMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( TRANS.EQ.ILATRANS( 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Set SAFE1 essentially to be the underflow threshold times the
*     number of additions in each row.
*
      SAFE1 = SLAMCH( 'Safe minimum' )
      SAFE1 = (N+1)*SAFE1
*
*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y).
*
*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
*     the inexact flag.  Still doesn't help change the iteration order
*     to per-column.
*
      IY = KY
      IF ( INCX.EQ.1 ) THEN
         IF( TRANS.EQ.ILATRANS( 'N' ) )THEN
            DO I = 1, LENY
               IF ( BETA .EQ. 0.0 ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               ELSE IF ( Y( IY ) .EQ. 0.0 ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. 0.0 ) THEN
                  DO J = 1, LENX
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO ) Y( IY ) =
     $              Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         ELSE
            DO I = 1, LENY
               IF ( BETA .EQ. 0.0 ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               ELSE IF ( Y( IY ) .EQ. 0.0 ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. 0.0 ) THEN
                  DO J = 1, LENX
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( J ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( J ) )*TEMP
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO ) Y( IY ) =
     $              Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         END IF
      ELSE
         IF( TRANS.EQ.ILATRANS( 'N' ) )THEN
            DO I = 1, LENY
               IF ( BETA .EQ. 0.0 ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               ELSE IF ( Y( IY ) .EQ. 0.0 ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. 0.0 ) THEN
                  JX = KX
                  DO J = 1, LENX
                     TEMP = CABS1( A( I, J ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( JX ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO ) Y( IY ) =
     $              Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         ELSE
            DO I = 1, LENY
               IF ( BETA .EQ. 0.0 ) THEN
                  SYMB_ZERO = .TRUE.
                  Y( IY ) = 0.0
               ELSE IF ( Y( IY ) .EQ. 0.0 ) THEN
                  SYMB_ZERO = .TRUE.
               ELSE
                  SYMB_ZERO = .FALSE.
                  Y( IY ) = BETA * ABS( Y( IY ) )
               END IF
               IF ( ALPHA .NE. 0.0 ) THEN
                  JX = KX
                  DO J = 1, LENX
                     TEMP = CABS1( A( J, I ) )
                     SYMB_ZERO = SYMB_ZERO .AND.
     $                    ( X( JX ) .EQ. ZERO .OR. TEMP .EQ. ZERO )

                     Y( IY ) = Y( IY ) + ALPHA*CABS1( X( JX ) )*TEMP
                     JX = JX + INCX
                  END DO
               END IF

               IF ( .NOT.SYMB_ZERO ) Y( IY ) =
     $              Y( IY ) + SIGN( SAFE1, Y( IY ) )

               IY = IY + INCY
            END DO
         END IF

      END IF
*
      RETURN
*
*     End of CLA_GEAMV
*
      END
