      subroutine orth_random ( n, seed, a )

c*********************************************************************72
c
cc ORTH_RANDOM returns the ORTH_RANDOM matrix.
c
c  Discussion:
c
c    The matrix is a random orthogonal matrix.
c
c  Properties:
c
c    The inverse of A is equal to A'.
c    A is orthogonal: A * A' = A' * A = I.
c    Because A is orthogonal, it is normal: A' * A = A * A'.
c    Columns and rows of A have unit Euclidean norm.
c    Distinct pairs of columns of A are orthogonal.
c    Distinct pairs of rows of A are orthogonal.
c    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
c    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
c    det ( A ) = +1 or -1.
c    A is unimodular.
c    All the eigenvalues of A have modulus 1.
c    All singular values of A are 1.
c    All entries of A are between -1 and 1.
c
c  Discussion:
c
c    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
c    National Academy of Sciences of Belarus, for convincingly
c    pointing out the severe deficiencies of an earlier version of
c    this routine.
c
c    Essentially, the computation involves saving the Q factor of the
c    QR factorization of a matrix whose entries are normally distributed.
c    However, it is only necessary to generate this matrix a column at
c    a time, since it can be shown that when it comes time to annihilate
c    the subdiagonal elements of column K, these (transformed) elements of
c    column K are still normally distributed random values.  Hence, there
c    is no need to generate them at the beginning of the process and
c    transform them K-1 times.
c
c    For computational efficiency, the individual Householder transformations
c    could be saved, as recommended in the reference, instead of being
c    accumulated into an explicit matrix format.
c
c  Licensing:
c
c    This code is distributed under the MIT license.
c
c  Modified:
c
c    By Sourangshu Ghosh
c
c  Author:
c
c    Sourangshu Ghosh
c
c  Reference:
c
c    Pete Stewart,
c    Efficient Generation of Random Orthogonal Matrices With an Application
c    to Condition Estimators,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 3, June 1980, pages 403-409.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision A(N,N), the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j
      double precision r8_normal_01
      integer seed
      double precision v(n)
      double precision x(n)
c
c  Start with A = the identity matrix.
c
      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            a(i,j) = 1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do
c
c  Now behave as though we were computing the QR factorization of
c  some other random matrix.  Generate the N elements of the first column,
c  compute the Householder matrix H1 that annihilates the subdiagonal elements,
c  and set A := A * H1' = A * H.
c
c  On the second step, generate the lower N-1 elements of the second column,
c  compute the Householder matrix H2 that annihilates them,
c  and set A := A * H2' = A * H2 = H1 * H2.
c
c  On the N-1 step, generate the lower 2 elements of column N-1,
c  compute the Householder matrix HN-1 that annihilates them, and
c  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
c  This is our random orthogonal matrix.
c
      do j = 1, n - 1
c
c  Set the vector that represents the J-th column to be annihilated.
c
        do i = 1, j - 1
          x(i) = 0.0D+00
        end do

        do i = j, n
          x(i) = r8_normal_01 ( seed )
        end do
c
c  Compute the vector V that defines a Householder transformation matrix
c  H(V) that annihilates the subdiagonal elements of X.
c
        call r8vec_house_column ( n, x, j, v )
c
c  Postmultiply the matrix A by H'(V) = H(V).
c
        call r8mat_house_axh ( n, a, v, a )

      end do

      return
      end
      subroutine pds_random ( n, seed, a )

c*********************************************************************72
c
cc PDS_RANDOM returns the PDS_RANDOM matrix.
c
c  Discussion:
c
c    The matrix is a "random" positive definite symmetric matrix.
c
c    The matrix returned will have eigenvalues in the range [0,1].
c
c  Properties:
c
c    A is symmetric: A' = A.
c
c    A is positive definite: 0 .lt. x'*A*x for nonzero x.
c
c    The eigenvalues of A will be real.
c
c  Licensing:
c
c    This code is distributed under the MIT license.
c
c  Modified:
c
c    BY Sourangshu Ghosh
c
c  Author:
c
c    Sourangshu Ghosh
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A(N,N), the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j
      integer k
      double precision lambda(n)
      double precision q(n,n)
      integer seed
c
c  Get a random set of eigenvalues.
c
      call r8vec_uniform_01 ( n, seed, lambda )
c
c  Get a random orthogonal matrix Q.
c
      call orth_random ( n, seed, q )
c
c  Set A = Q * Lambda * Q'.
c
      do i = 1, n
        do j = 1, n
          a(i,j) = 0.0D+00
          do k = 1, n
            a(i,j) = a(i,j) + q(i,k) * lambda(k) * q(j,k)
          end do
        end do
      end do

      return
      end
      function r8_normal_01 ( seed )

c*********************************************************************72
c
cc R8_NORMAL_01 returns a unit pseudonormal R8.
c
c  Discussion:
c
c    Because this routine uses the Box Muller method, it requires pairs
c    of uniform random values to generate a pair of normal random values.
c    This means that on every other call, the code can use the second
c    value that it calculated.
c
c    However, if the user has changed the SEED value between calls,
c    the routine automatically resets itself and discards the saved data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    Sourangshu Ghosh
c
c  Author:
c
c    Sourangshu Ghosh
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision R8_NORMAL_01, a sample of the standard normal PDF.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_normal_01
      double precision r8_uniform_01
      integer seed
      integer seed1
      integer seed2
      integer seed3
      integer used
      double precision v1
      double precision v2

      save seed1
      save seed2
      save seed3
      save used
      save v2

      data seed2 / 0 /
      data used / 0 /
      data v2 / 0.0D+00 /
c
c  If USED is odd, but the input SEED does not match
c  the output SEED on the previous call, then the user has changed
c  the seed.  Wipe out internal memory.
c
      if ( mod ( used, 2 ) .eq. 1 ) then

        if ( seed .ne. seed2 ) then
          used = 0
          seed1 = 0
          seed2 = 0
          seed3 = 0
          v2 = 0.0D+00
        end if

      end if
c
c  If USED is even, generate two uniforms, create two normals,
c  return the first normal and its corresponding seed.
c
      if ( mod ( used, 2 ) .eq. 0 ) then

        seed1 = seed

        r1 = r8_uniform_01 ( seed )

        if ( r1 .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop 1
        end if

        seed2 = seed

        r2 = r8_uniform_01 ( seed )

        seed3 = seed

        v1 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
        v2 = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )

        r8_normal_01 = v1
        seed = seed2
c
c  If USED is odd (and the input SEED matched the output value from
c  the previous call), return the second normal and its corresponding seed.
c
      else

        r8_normal_01 = v2
        seed = seed3

      end if

      used = used + 1

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r8_uniform_01 = seed / ( 2^31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    By Sourangshu Ghosh
c
c  Author:
c
c    Sourangshu Ghosh
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop 1
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r83_cg ( n, a, b, x )

c*********************************************************************72
c
cc R83_CG uses the conjugate gradient method on an R83 system.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    The matrix A must be a positive definite symmetric band matrix.
c
c    The method is designed to reach the solution after N computational
c    steps.  However, roundoff may introduce unacceptably large errors for
c    some problems.  In such a case, calling the routine again, using
c    the computed solution as the new starting estimate, should improve
c    the results.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    By Sourangshu Ghosh
c
c  Author:
c
c    Sourangshu Ghosh
c
c  Reference:
c
c    Frank Beckman,
c    The Solution of Linear Equations by the Conjugate Gradient Method,
c    in Mathematical Methods for Digital Computers,
c    edited by John Ralston, Herbert Wilf,
c    Wiley, 1967,
c    ISBN: 0471706892,
c    LC: QA76.5.R3.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(3,N), the matrix.
c
c    Input, double precision B(N), the right hand side vector.
c
c    Input/output, double precision X(N).
c    On input, an estimate for the solution, which may be 0.
c    On output, the approximate solution vector.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      integer i
      integer it
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      double precision r8vec_dot_product
      double precision rap
      double precision x(n)
c
c  Initialize
c    AP = A * x,
c    R  = b - A * x,
c    P  = b - A * x.
c
      call r83_mv ( n, n, a, x, ap )

      do i = 1, n
        r(i) = b(i) - ap(i)
      end do

      do i = 1, n
        p(i) = b(i) - ap(i)
      end do
c
c  Do the N steps of the conjugate gradient method.
c
      do it = 1, n
c
c  Compute the matrix*vector product AP=A*P.
c
        call r83_mv ( n, n, a, p, ap )
c
c  Compute the dot products
c    PAP = P*AP,
c    PR  = P*R
c  Set
c    ALPHA = PR / PAP.
c
        pap = r8vec_dot_product ( n, p, ap )
        pr = r8vec_dot_product ( n, p, r )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
c
c  Set
c    X = X + ALPHA * P
c    R = R - ALPHA * AP.
c
        do i = 1, n
          x(i) = x(i) + alpha * p(i)
        end do

        do i = 1, n
          r(i) = r(i) - alpha * ap(i)
        end do
c
c  Compute the vector dot product
c    RAP = R*AP
c  Set
c    BETA = - RAP / PAP.
c
        rap = r8vec_dot_product ( n, r, ap )

        beta = - rap / pap
c
c  Update the perturbation vector
c    P = R + BETA * P.
c
        do i = 1, n
          p(i) = r(i) + beta * p(i)
        end do

      end do

      return
      end
      subroutine r83_dif2 ( m, n, a )

c*********************************************************************72
c
cc R83_DIF2 returns the DIF2 matrix in R83 format.
c
c  Example:
c
c    N = 5
c
c    2 -1  .  .  .
c   -1  2 -1  .  .
c    . -1  2 -1  .
c    .  . -1  2 -1
c    .  .  . -1  2
c
c  Properties:
c
c    A is banded, with bandwidth 3.
c
c    A is tridiagonal.
c
c    Because A is tridiagonal, it has property A (bipartite).
c
c    A is a special case of the TRIS or tridiagonal scalar matrix.
c
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c
c    A is Toeplitz: constant along diagonals.
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    A is positive definite.
c
c    A is an M matrix.
c
c    A is weakly diagonally dominant, but not strictly diagonally dominant.
c
c    A has an LU factorization A = L * U, without pivoting.
c
c      The matrix L is lower bidiagonal with subdiagonal elements:
c
c        L(I+1,I) = -I/(I+1)
c
c      The matrix U is upper bidiagonal, with diagonal elements
c
c        U(I,I) = (I+1)/I
c
c      and superdiagonal elements which are all -1.
c
c    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
c
c      L(I,I) =    sqrt ( (I+1) / I )
c      L(I,I-1) = -sqrt ( (I-1) / I )
c
c    The eigenvalues are
c
c      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
c                = 4 SIN^2(I*PI/(2*N+2))
c
c    The corresponding eigenvector X(I) has entries
c
c       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
c
c    Simple linear systems:
c
c      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
c
c      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
c
c    det ( A ) = N + 1.
c
c    The value of the determinant can be seen by induction,
c    and expanding the determinant across the first row:
c
c      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
c                = 2 * N - (N-1)
c                = N + 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c   Sourangshu Ghosh
c
c  Author:
c
c    Sourangshu Ghosh
c
c  Reference:
c
c    Robert Gregory, David Karney,
c    A Collection of Matrices for Testing Computational Algorithms,
c    Wiley, 1969,
c    ISBN: 0882756494,
c    LC: QA263.68
c
c    Morris Newman, John Todd,
c    Example A8,
c    The evaluation of matrix inversion programs,
c    Journal of the Society for Industrial and Applied Mathematics,
c    Volume 6, Number 4, pages 466-476, 1958.
c
c    John Todd,
c    Basic Numerical Mathematics,
c    Volume 2: Numerical Algebra,
c    Birkhauser, 1980,
c    ISBN: 0817608117,
c    LC: QA297.T58.
c
c    Joan Westlake,
c    A Handbook of Numerical Matrix Inversion and Solution of 
c    Linear Equations,
c    John Wiley, 1968,
c    ISBN13: 978-0471936756,
c    LC: QA263.W47.
c
c  Parameters:
c
c    Input, integer M, N, the order of the matrix.
c
c    Output, double precision A(3,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(3,n)
      integer i
      integer j
      integer mn

      do j = 1, n
        do i = 1, 3
          a(i,j) = 0.0D+00
        end do
      end do

      mn = min ( m, n )

      do j = 1, mn
        a(2,j)   = +2.0D+00
      end do

      do j = 2, mn
        a(1,j)   = -1.0D+00
      end do

      if ( m .le. n ) then
        do j = 1, mn - 1
          a(3,j) = -1.0D+00
        end do
      else if ( n .lt. m ) then
        do j = 1, mn
          a(3,j) = -1.0D+00
        end do
      end if
      
      return
      end
      subroutine r83_mv ( m, n, a, x, b )

c*********************************************************************72
c
cc R83_MV multiplies an R83 matrix times an R8VEC.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(3,N), the R83 matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer m
      integer n

      double precision a(3,n)
      double precision b(m)
      integer i
      integer j
      integer mn
      double precision x(n)

      do i = 1, m
        b(i) = 0.0D+00
      end do

      mn = min ( m, n )

      if ( n .eq. 1 ) then
        b(1) = a(2,1) * x(1)
        if ( 1 .lt. m ) then
          b(2) = a(3,1) * x(1)
        end if
        return
      end if

      b(1) = a(2,1)      * x(1) 
     &     + a(1,2)      * x(2)

      do j = 2, mn - 1
        b(j) = a(3,j-1) * x(j-1) 
     &       + a(2,j) * x(j) 
     &       + a(1,j+1)   * x(j+1)
      end do

      b(mn)     = a(3,mn-1)   * x(mn-1) 
     &          + a(2,mn)     * x(mn)

      if ( n .lt. m ) then
        b(mn+1) = b(mn+1) + a(3,mn) * x(mn)
      end if

      if ( m .lt. n ) then
        b(mn) = b(mn) + a(1,mn+1) * x(mn+1)
      end if

      return
      end
      subroutine r83_resid ( m, n, a, x, b, r )

c*********************************************************************72
c
cc R83_RESID computes the residual R = B-A*X for R83 matrices.
c
c  Licensing:
c
c    This code is distributed under the MIT license. 
c
c  Modified:
c
c    02 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(3,N), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Input, double precision B(M), the desired result A * x.
c
c    Output, double precision R(M), the residual R = B - A * X.
c
      implicit none

      integer m
      integer n

      double precision a(3,n)
      double precision b(m)
      integer i
      double precision r(m)
      double precision x(n)

      call r83_mv ( m, n, a, x, r )

      do i = 1, m
        r(i) = b(i) - r(i)
      end do

      return
      end
      subroutine r83s_cg ( n, a, b, x )

c*********************************************************************72
c
cc R83S_CG uses the conjugate gradient method on an R83S system.
c
c  Discussion:
c
c    The R83S storage format is used for a tridiagonal scalar matrix.
c    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
c    values that occur on every row.
c
c    The matrix A must be a positive definite symmetric band matrix.
c
c    The method is designed to reach the solution after N computational
c    steps.  However, roundoff may introduce unacceptably large errors for
c    some problems.  In such a case, calling the routine again, using
c    the computed solution as the new starting estimate, should improve
c    the results.
c
c  Example:
c
c    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
c    be interpreted:
c
c      A2  A3   0   0   0
c      A1  A2  A3   0   0
c       0  A1  A2  A3   0 
c       0   0  A1  A2  A3
c       0   0   0  A1  A2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 July 2014
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Frank Beckman,
c    The Solution of Linear Equations by the Conjugate Gradient Method,
c    in Mathematical Methods for Digital Computers,
c    edited by John Ralston, Herbert Wilf,
c    Wiley, 1967,
c    ISBN: 0471706892,
c    LC: QA76.5.R3.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(3), the matrix.
c
c    Input, double precision B(N), the right hand side vector.
c
c    Input/output, double precision X(N).
c    On input, an estimate for the solution, which may be 0.
c    On output, the approximate solution vector.
c
      implicit none

      integer n

      double precision a(3)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      integer i
      integer it
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      double precision r8vec_dot_product
      double precision rap
      double precision x(n)
c
c  Initialize
c    AP = A * x,
c    R  = b - A * x,
c    P  = b - A * x.
c
      call r83s_mv ( n, n, a, x, ap )

      do i = 1, n
        r(i) = b(i) - ap(i)
      end do

      do i = 1, n
        p(i) = b(i) - ap(i)
      end do
c
c  Do the N steps of the conjugate gradient method.
c
      do it = 1, n
c
c  Compute the matrix*vector product AP=A*P.
c
        call r83s_mv ( n, n, a, p, ap )
c
c  Compute the dot products
c    PAP = P*AP,
c    PR  = P*R
c  Set
c    ALPHA = PR / PAP.
c
        pap = r8vec_dot_product ( n, p, ap )
        pr = r8vec_dot_product ( n, p, r )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
c
c  Set
c    X = X + ALPHA * P
c    R = R - ALPHA * AP.
c
        do i = 1, n
          x(i) = x(i) + alpha * p(i)
        end do

        do i = 1, n
          r(i) = r(i) - alpha * ap(i)
        end do
c
c  Compute the vector dot product
c    RAP = R*AP
c  Set
c    BETA = - RAP / PAP.
c
        rap = r8vec_dot_product ( n, r, ap )

        beta = - rap / pap
c
c  Update the perturbation vector
c    P = R + BETA * P.
c
        do i = 1, n
          p(i) = r(i) + beta * p(i)
        end do

      end do

      return
      end
      subroutine r83s_dif2 ( m, n, a )

c*********************************************************************72
c
cc R83S_DIF2 returns the DIF2 matrix in R83S format.
c
c  Example:
c
c    N = 5
c
c    2 -1  .  .  .
c   -1  2 -1  .  .
c    . -1  2 -1  .
c    .  . -1  2 -1
c    .  .  . -1  2
c
c  Properties:
c
c    A is banded, with bandwidth 3.
c    A is tridiagonal.
c    Because A is tridiagonal, it has property A (bipartite).
c    A is a special case of the TRIS or tridiagonal scalar matrix.
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c    A is Toeplitz: constant along diagonals.
c    A is symmetric: A' = A.
c    Because A is symmetric, it is normal.
c    Because A is normal, it is diagonalizable.
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c    A is positive definite.
c    A is an M matrix.
c    A is weakly diagonally dominant, but not strictly diagonally dominant.
c    A has an LU factorization A = L * U, without pivoting.
c      The matrix L is lower bidiagonal with subdiagonal elements:
c        L(I+1,I) = -I/(I+1)
c      The matrix U is upper bidiagonal, with diagonal elements
c        U(I,I) = (I+1)/I
c      and superdiagonal elements which are all -1.
c    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
c      L(I,I) =    sqrt ( (I+1) / I )
c      L(I,I-1) = -sqrt ( (I-1) / I )
c    The eigenvalues are
c      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
c                = 4 SIN^2(I*PI/(2*N+2))
c    The corresponding eigenvector X(I) has entries
c       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
c    Simple linear systems:
c      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
c      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
c    det ( A ) = N + 1.
c    The value of the determinant can be seen by induction,
c    and expanding the determinant across the first row:
c      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
c                = 2 * N - (N-1)
c                = N + 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 July 2014
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Gregory, David Karney,
c    A Collection of Matrices for Testing Computational Algorithms,
c    Wiley, 1969,
c    ISBN: 0882756494,
c    LC: QA263.68
c
c    Morris Newman, John Todd,
c    Example A8,
c    The evaluation of matrix inversion programs,
c    Journal of the Society for Industrial and Applied Mathematics,
c    Volume 6, Number 4, pages 466-476, 1958.
c
c    John Todd,
c    Basic Numerical Mathematics,
c    Volume 2: Numerical Algebra,
c    Birkhauser, 1980,
c    ISBN: 0817608117,
c    LC: QA297.T58.
c
c    Joan Westlake,
c    A Handbook of Numerical Matrix Inversion and Solution of 
c    Linear Equations,
c    John Wiley, 1968,
c    ISBN13: 978-0471936756,
c    LC: QA263.W47.
c
c  Parameters:
c
c    Input, integer M, N, the order of the matrix.
c
c    Output, double precision A(3), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(3)

      a(1) = -1.0D+00
      a(2) = +2.0D+00
      a(3) = -1.0D+00

      return
      end
      subroutine r83s_mv ( m, n, a, x, b )

c*********************************************************************72
c
cc R83S_MV multiplies an R83S matrix times an R8VEC.
c
c  Discussion:
c
c    The R83S storage format is used for a tridiagonal scalar matrix.
c    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
c    values that occur on every row.
c
c  Example:
c
c    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
c    be interpreted:
c
c      A2  A3   0   0   0
c      A1  A2  A3   0   0
c       0  A1  A2  A3   0 
c       0   0  A1  A2  A3
c       0   0   0  A1  A2
c
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 July 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(3), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer m
      integer n

      double precision a(3)
      double precision b(m)
      integer i
      double precision x(n)

      do i = 1, m
        b(i) = 0.0D+00
      end do

      do i = 2, n
        b(i) = b(i) + a(1) * x(i-1)
      end do

      do i = 1, n
        b(i) = b(i) + a(2) * x(i)
      end do

      do i = 1, n - 1
        b(i) = b(i) + a(3) * x(i+1)
      end do

      return
      end
      subroutine r83s_resid ( m, n, a, x, b, r )

c*********************************************************************72
c
cc R83S_RESID computes the residual R = B-A*X for R83S matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 July 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(3), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Input, double precision B(M), the desired result A * x.
c
c    Output, double precision R(M), the residual R = B - A * X.
c
      implicit none

      integer m
      integer n

      double precision a(3)
      double precision b(m)
      integer i
      double precision r(m)
      double precision x(n)

      call r83s_mv ( m, n, a, x, r )

      do i = 1, m
        r(i) = b(i) - r(i)
      end do

      return
      end
      subroutine r83t_cg ( n, a, b, x )

c*********************************************************************72
c
cc R83T_CG uses the conjugate gradient method on an R83T system.
c
c  Discussion:
c
c    The R83T storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1:N-1,3), the diagonal in
c    entries (1:N,2), and the subdiagonal in (2:N,1).  Thus, the
c    original matrix is "collapsed" horizontally into the array.
c
c    The matrix A must be a positive definite symmetric band matrix.
c
c    The method is designed to reach the solution after N computational
c    steps.  However, roundoff may introduce unacceptably large errors for
c    some problems.  In such a case, calling the routine again, using
c    the computed solution as the new starting estimate, should improve
c    the results.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Frank Beckman,
c    The Solution of Linear Equations by the Conjugate Gradient Method,
c    in Mathematical Methods for Digital Computers,
c    edited by John Ralston, Herbert Wilf,
c    Wiley, 1967,
c    ISBN: 0471706892,
c    LC: QA76.5.R3.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(N,3), the matrix.
c
c    Input, double precision B(N), the right hand side vector.
c
c    Input/output, double precision X(N).
c    On input, an estimate for the solution, which may be 0.
c    On output, the approximate solution vector.
c
      implicit none

      integer n

      double precision a(n,3)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      integer i
      integer it
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      double precision r8vec_dot_product
      double precision rap
      double precision x(n)
c
c  Initialize
c    AP = A * x,
c    R  = b - A * x,
c    P  = b - A * x.
c
      call r83t_mv ( n, n, a, x, ap )

      do i = 1, n
        r(i) = b(i) - ap(i)
      end do

      do i = 1, n
        p(i) = b(i) - ap(i)
      end do
c
c  Do the N steps of the conjugate gradient method.
c
      do it = 1, n
c
c  Compute the matrix*vector product AP=A*P.
c
        call r83t_mv ( n, n, a, p, ap )
c
c  Compute the dot products
c    PAP = P*AP,
c    PR  = P*R
c  Set
c    ALPHA = PR / PAP.
c
        pap = r8vec_dot_product ( n, p, ap )
        pr = r8vec_dot_product ( n, p, r )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
c
c  Set
c    X = X + ALPHA * P
c    R = R - ALPHA * AP.
c
        do i = 1, n
          x(i) = x(i) + alpha * p(i)
        end do

        do i = 1, n
          r(i) = r(i) - alpha * ap(i)
        end do
c
c  Compute the vector dot product
c    RAP = R*AP
c  Set
c    BETA = - RAP / PAP.
c
        rap = r8vec_dot_product ( n, r, ap )

        beta = - rap / pap
c
c  Update the perturbation vector
c    P = R + BETA * P.
c
        do i = 1, n
          p(i) = r(i) + beta * p(i)
        end do

      end do

      return
      end
      subroutine r83t_dif2 ( m, n, a )

c*********************************************************************72
c
cc R83T_DIF2 returns the DIF2 matrix in R83T format.
c
c  Example:
c
c    N = 5
c
c    2 -1  .  .  .
c   -1  2 -1  .  .
c    . -1  2 -1  .
c    .  . -1  2 -1
c    .  .  . -1  2
c
c  Properties:
c
c    A is banded, with bandwidth 3.
c
c    A is tridiagonal.
c
c    Because A is tridiagonal, it has property A (bipartite).
c
c    A is a special case of the TRIS or tridiagonal scalar matrix.
c
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c
c    A is Toeplitz: constant along diagonals.
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    A is positive definite.
c
c    A is an M matrix.
c
c    A is weakly diagonally dominant, but not strictly diagonally dominant.
c
c    A has an LU factorization A = L * U, without pivoting.
c
c      The matrix L is lower bidiagonal with subdiagonal elements:
c
c        L(I+1,I) = -I/(I+1)
c
c      The matrix U is upper bidiagonal, with diagonal elements
c
c        U(I,I) = (I+1)/I
c
c      and superdiagonal elements which are all -1.
c
c    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
c
c      L(I,I) =    sqrt ( (I+1) / I )
c      L(I,I-1) = -sqrt ( (I-1) / I )
c
c    The eigenvalues are
c
c      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
c                = 4 SIN^2(I*PI/(2*N+2))
c
c    The corresponding eigenvector X(I) has entries
c
c       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
c
c    Simple linear systems:
c
c      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
c
c      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
c
c    det ( A ) = N + 1.
c
c    The value of the determinant can be seen by induction,
c    and expanding the determinant across the first row:
c
c      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
c                = 2 * N - (N-1)
c                = N + 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Gregory, David Karney,
c    A Collection of Matrices for Testing Computational Algorithms,
c    Wiley, 1969,
c    ISBN: 0882756494,
c    LC: QA263.68
c
c    Morris Newman, John Todd,
c    Example A8,
c    The evaluation of matrix inversion programs,
c    Journal of the Society for Industrial and Applied Mathematics,
c    Volume 6, Number 4, pages 466-476, 1958.
c
c    John Todd,
c    Basic Numerical Mathematics,
c    Volume 2: Numerical Algebra,
c    Birkhauser, 1980,
c    ISBN: 0817608117,
c    LC: QA297.T58.
c
c    Joan Westlake,
c    A Handbook of Numerical Matrix Inversion and Solution of 
c    Linear Equations,
c    John Wiley, 1968,
c    ISBN13: 978-0471936756,
c    LC: QA263.W47.
c
c  Parameters:
c
c    Input, integer M, N, the order of the matrix.
c
c    Output, double precision A(M,3), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,3)
      integer i
      integer j
      integer mn

      do j = 1, 3
        do i = 1, m
          a(i,j) = 0.0D+00
        end do
      end do

      mn = min ( m, n )

      do i = 2, mn
        a(i,1) = -1.0D+00
      end do

      do i = 1, mn
        a(i,2) = +2.0D+00
      end do

      do i = 1, mn - 1
        a(i,3) = -1.0D+00
      end do

      if ( m .lt. n ) then
        a(mn,3) = -1.0D+00
      else if ( n .lt. m ) then
        a(mn+1,1) = -1.0D+00
      end if
      
      return
      end
      subroutine r83t_mv ( m, n, a, x, b )

c*********************************************************************72
c
cc R83T_MV multiplies an R83T matrix times an R8VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,3), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer m
      integer n

      double precision a(m,3)
      double precision b(m)
      integer i
      integer mn
      double precision x(n)

      do i = 1, m
        b(i) = 0.0D+00
      end do

      mn = min ( m, n )

      if ( n .eq. 1 ) then
        b(1) = a(1,2) * x(1)
        if ( 1 .lt. m ) then
          b(2) = a(2,1) * x(1)
        end if
        return
      end if

      b(1) = a(1,2)      * x(1) 
     &     + a(1,3)      * x(2)

      do i = 2, mn - 1
        b(i) = a(i,1) * x(i-1) 
     &       + a(i,2) * x(i) 
     &       + a(i,3) * x(i+1)
      end do

      b(mn)     = a(mn,1) * x(mn-1) 
     &          + a(mn,2) * x(mn)

      if ( n .lt. m ) then
        b(mn+1) = b(mn+1) + a(mn+1,1) * x(mn)
      else if ( m .lt. n ) then
        b(mn) = b(mn) + a(mn,3) * x(mn+1)
      end if

      return
      end
      subroutine r83t_resid ( m, n, a, x, b, r )

c*********************************************************************72
c
cc R83T_RESID computes the residual R = B-A*X for R83T matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,3), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Input, double precision B(M), the desired result A * x.
c
c    Output, double precision R(M), the residual R = B - A * X.
c
      implicit none

      integer m
      integer n

      double precision a(m,3)
      double precision b(m)
      integer i
      double precision r(m)
      double precision x(n)

      call r83t_mv ( m, n, a, x, r )

      do i = 1, m
        r(i) = b(i) - r(i)
      end do

      return
      end
      subroutine r8ge_cg ( n, a, b, x )

c*********************************************************************72
c
cc R8GE_CG uses the conjugate gradient method on an R8GE system.
c
c  Discussion:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c    The matrix A must be a positive definite symmetric band matrix.
c
c    The method is designed to reach the solution after N computational
c    steps.  However, roundoff may introduce unacceptably large errors for
c    some problems.  In such a case, calling the routine again, using
c    the computed solution as the new starting estimate, should improve
c    the results.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Frank Beckman,
c    The Solution of Linear Equations by the Conjugate Gradient Method,
c    in Mathematical Methods for Digital Computers,
c    edited by John Ralston, Herbert Wilf,
c    Wiley, 1967,
c    ISBN: 0471706892,
c    LC: QA76.5.R3.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(N,N), the matrix.
c
c    Input, double precision B(N), the right hand side vector.
c
c    Input/output, double precision X(N).
c    On input, an estimate for the solution, which may be 0.
c    On output, the approximate solution vector.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      integer i
      integer it
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      double precision r8vec_dot_product
      double precision rap
      double precision x(n)
c
c  Initialize
c    AP = A * x,
c    R  = b - A * x,
c    P  = b - A * x.
c
      call r8ge_mv ( n, n, a, x, ap )

      do i = 1, n 
        r(i) = b(i) - ap(i)
      end do

      do i = 1, n
        p(i) = b(i) - ap(i)
      end do
c
c  Do the N steps of the conjugate gradient method.
c
      do it = 1, n
c
c  Compute the matrix*vector product AP=A*P.
c
        call r8ge_mv ( n, n, a, p, ap )
c
c  Compute the dot products
c    PAP = P*AP,
c    PR  = P*R
c  Set
c    ALPHA = PR / PAP.
c
        pap = r8vec_dot_product ( n, p, ap )
        pr = r8vec_dot_product ( n, p, r )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
c
c  Set
c    X = X + ALPHA * P
c    R = R - ALPHA * AP.
c
        do i = 1, n
          x(i) = x(i) + alpha * p(i)
        end do

        do i = 1, n
          r(i) = r(i) - alpha * ap(i)
        end do
c
c  Compute the vector dot product
c    RAP = R*AP
c  Set
c    BETA = - RAP / PAP.
c
        rap = r8vec_dot_product ( n, r, ap )

        beta = - rap / pap
c
c  Update the perturbation vector
c    P = R + BETA * P.
c
        do i = 1, n
          p(i) = r(i) + beta * p(i)
        end do

      end do

      return
      end
      subroutine r8ge_dif2 ( m, n, a )

c*********************************************************************72
c
cc R8GE_DIF2 returns the DIF2 matrix in R8GE format.
c
c  Example:
c
c    N = 5
c
c    2 -1  .  .  .
c   -1  2 -1  .  .
c    . -1  2 -1  .
c    .  . -1  2 -1
c    .  .  . -1  2
c
c  Properties:
c
c    A is banded, with bandwidth 3.
c
c    A is tridiagonal.
c
c    Because A is tridiagonal, it has property A (bipartite).
c
c    A is a special case of the TRIS or tridiagonal scalar matrix.
c
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c
c    A is Toeplitz: constant along diagonals.
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    A is positive definite.
c
c    A is an M matrix.
c
c    A is weakly diagonally dominant, but not strictly diagonally dominant.
c
c    A has an LU factorization A = L * U, without pivoting.
c
c      The matrix L is lower bidiagonal with subdiagonal elements:
c
c        L(I+1,I) = -I/(I+1)
c
c      The matrix U is upper bidiagonal, with diagonal elements
c
c        U(I,I) = (I+1)/I
c
c      and superdiagonal elements which are all -1.
c
c    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
c
c      L(I,I) =    sqrt ( (I+1) / I )
c      L(I,I-1) = -sqrt ( (I-1) / I )
c
c    The eigenvalues are
c
c      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
c                = 4 SIN^2(I*PI/(2*N+2))
c
c    The corresponding eigenvector X(I) has entries
c
c       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
c
c    Simple linear systems:
c
c      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
c
c      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
c
c    det ( A ) = N + 1.
c
c    The value of the determinant can be seen by induction,
c    and expanding the determinant across the first row:
c
c      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
c                = 2 * N - (N-1)
c                = N + 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 July 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Gregory, David Karney,
c    A Collection of Matrices for Testing Computational Algorithms,
c    Wiley, 1969,
c    ISBN: 0882756494,
c    LC: QA263.68
c
c    Morris Newman, John Todd,
c    Example A8,
c    The evaluation of matrix inversion programs,
c    Journal of the Society for Industrial and Applied Mathematics,
c    Volume 6, Number 4, pages 466-476, 1958.
c
c    John Todd,
c    Basic Numerical Mathematics,
c    Volume 2: Numerical Algebra,
c    Birkhauser, 1980,
c    ISBN: 0817608117,
c    LC: QA297.T58.
c
c    Joan Westlake,
c    A Handbook of Numerical Matrix Inversion and Solution of 
c    Linear Equations,
c    John Wiley, 1968,
c    ISBN13: 978-0471936756,
c    LC: QA263.W47.
c
c  Parameters:
c
c    Input, integer M, N, the order of the matrix.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m

          if ( j .eq. i - 1 ) then
            a(i,j) = -1.0D+00
          else if ( j .eq. i ) then
            a(i,j) = 2.0D+00
          else if ( j .eq. i + 1 ) then
            a(i,j) = -1.0D+00
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do

      return
      end
      subroutine r8ge_mv ( m, n, a, x, b )

c*********************************************************************72
c
cc R8GE_MV multiplies an R8GE matrix by an R8VEC.
c
c  Discussion:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c    R8GE storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8GE matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(m)
      double precision x(n)

      call r8mat_mv ( m, n, a, x, b )

      return
      end
      subroutine r8ge_resid ( m, n, a, x, b, r )

c*********************************************************************72
c
cc R8GE_RESID computes the residual R = B-A*X for R8GE matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Input, double precision B(M), the desired result A * x.
c
c    Output, double precision R(M), the residual R = B - A * X.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(m)
      integer i
      double precision r(m)
      double precision x(n)

      call r8mat_mv ( m, n, a, x, r )

      do i = 1, m
        r(i) = b(i) - r(i)
      end do

      return
      end
      subroutine r8mat_house_axh ( n, a, v, ah )

c*********************************************************************72
c
cc R8MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
c
c  Discussion:
c
c    The Householder matrix H(V) is defined by
c
c      H(V) = I - 2 * v * v' / ( v' * v )
c
c    This routine is not particularly efficient.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of A.
c
c    Input, double precision A(N,N), the matrix.
c
c    Input, double precision V(N), a vector defining a Householder matrix.
c
c    Output, double precision AH(N,N), the product A*H.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision ah(n,n)
      double precision av(n)
      integer i
      integer j
      double precision v(n)
      double precision v_normsq

      v_normsq = 0.0D+00
      do i = 1, n
        v_normsq = v_normsq + v(i)**2
      end do
      
      do i = 1, n
        av(i) = 0.0D+00
        do j = 1, n
          av(i) = av(i) + a(i,j) * v(j)
        end do
      end do

      do i = 1, n
        do j = 1, n
          ah(i,j) = a(i,j)
        end do
      end do

      do i = 1, n
        do j = 1, n
          ah(i,j) = ah(i,j) - 2.0D+00 * av(i) * v(j)
        end do
      end do

      do i = 1, n
        do j = 1, n
          ah(i,j) = ah(i,j) / v_normsq
        end do
      end do

      return
      end
      subroutine r8mat_mv ( m, n, a, x, y )

c*********************************************************************72
c
cc R8MAT_MV multiplies a matrix times a vector.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c    In FORTRAN90, this operation can be more efficiently carried
c    out by the command
c
c      Y(1:M) = MATMUL ( A(1:M,1:N), X(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of the matrix.
c
c    Input, double precision A(M,N), the M by N matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision Y(M), the product A*X.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision x(n)
      double precision y(m)
      double precision y1(m)

      do i = 1, m
        y1(i) = 0.0D+00
        do j = 1, n
          y1(i) = y1(i) + a(i,j) * x(j)
        end do
      end do

      do i = 1, m
        y(i) = y1(i)
      end do

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 May 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r8pbu_cg ( n, mu, a, b, x )

c*********************************************************************72
c
cc R8PBU_CG uses the conjugate gradient method on an R8PBU system.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c    The matrix A must be a positive definite symmetric band matrix.
c
c    The method is designed to reach the solution after N computational
c    steps.  However, roundoff may introduce unacceptably large errors for
c    some problems.  In such a case, calling the routine again, using
c    the computed solution as the new starting estimate, should improve
c    the results.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Frank Beckman,
c    The Solution of Linear Equations by the Conjugate Gradient Method,
c    in Mathematical Methods for Digital Computers,
c    edited by John Ralston, Herbert Wilf,
c    Wiley, 1967,
c    ISBN: 0471706892,
c    LC: QA76.5.R3.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer MU, the number of superdiagonals.
c    MU must be at least 0, and no more than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBU matrix.
c
c    Input, double precision B(N), the right hand side vector.
c
c    Input/output, double precision X(N).
c    On input, an estimate for the solution, which may be 0.
c    On output, the approximate solution vector.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      integer i
      integer it
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      double precision r8vec_dot_product
      double precision rap
      double precision x(n)
c
c  Initialize
c    AP = A * x,
c    R  = b - A * x,
c    P  = b - A * x.
c
      call r8pbu_mv ( n, n, mu, a, x, ap )

      do i = 1, n 
        r(i) = b(i) - ap(i)
      end do

      do i = 1, n
        p(i) = b(i) - ap(i)
      end do
c
c  Do the N steps of the conjugate gradient method.
c
      do it = 1, n
c
c  Compute the matrix*vector product AP=A*P.
c
        call r8pbu_mv ( n, n, mu, a, p, ap )
c
c  Compute the dot products
c    PAP = P*AP,
c    PR  = P*R
c  Set
c    ALPHA = PR / PAP.
c
        pap = r8vec_dot_product ( n, p, ap )
        pr = r8vec_dot_product ( n, p, r )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
c
c  Set
c    X = X + ALPHA * P
c    R = R - ALPHA * AP.
c
        do i = 1, n
          x(i) = x(i) + alpha * p(i)
        end do

        do i = 1, n
          r(i) = r(i) - alpha * ap(i)
        end do
c
c  Compute the vector dot product
c    RAP = R*AP
c  Set
c    BETA = - RAP / PAP.
c
        rap = r8vec_dot_product ( n, r, ap )

        beta = - rap / pap
c
c  Update the perturbation vector
c    P = R + BETA * P.
c
        do i = 1, n
          p(i) = r(i) + beta * p(i)
        end do

      end do

      return
      end
      subroutine r8pbu_dif2 ( m, n, mu, a )

c*********************************************************************72
c
cc R8PBU_DIF2 returns the DIF2 matrix in R8PBU format.
c
c  Example:
c
c    N = 5
c
c    2 -1  .  .  .
c   -1  2 -1  .  .
c    . -1  2 -1  .
c    .  . -1  2 -1
c    .  .  . -1  2
c
c  Properties:
c
c    A is banded, with bandwidth 3.
c
c    A is tridiagonal.
c
c    Because A is tridiagonal, it has property A (bipartite).
c
c    A is a special case of the TRIS or tridiagonal scalar matrix.
c
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c
c    A is Toeplitz: constant along diagonals.
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    A is positive definite.
c
c    A is an M matrix.
c
c    A is weakly diagonally dominant, but not strictly diagonally dominant.
c
c    A has an LU factorization A = L * U, without pivoting.
c
c      The matrix L is lower bidiagonal with subdiagonal elements:
c
c        L(I+1,I) = -I/(I+1)
c
c      The matrix U is upper bidiagonal, with diagonal elements
c
c        U(I,I) = (I+1)/I
c
c      and superdiagonal elements which are all -1.
c
c    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
c
c      L(I,I) =    sqrt ( (I+1) / I )
c      L(I,I-1) = -sqrt ( (I-1) / I )
c
c    The eigenvalues are
c
c      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
c                = 4 SIN^2(I*PI/(2*N+2))
c
c    The corresponding eigenvector X(I) has entries
c
c       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
c
c    Simple linear systems:
c
c      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
c
c      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
c
c    det ( A ) = N + 1.
c
c    The value of the determinant can be seen by induction,
c    and expanding the determinant across the first row:
c
c      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
c                = 2 * N - (N-1)
c                = N + 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 July 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Gregory, David Karney,
c    A Collection of Matrices for Testing Computational Algorithms,
c    Wiley, 1969,
c    ISBN: 0882756494,
c    LC: QA263.68
c
c    Morris Newman, John Todd,
c    Example A8,
c    The evaluation of matrix inversion programs,
c    Journal of the Society for Industrial and Applied Mathematics,
c    Volume 6, Number 4, pages 466-476, 1958.
c
c    John Todd,
c    Basic Numerical Mathematics,
c    Volume 2: Numerical Algebra,
c    Birkhauser, 1980,
c    ISBN: 0817608117,
c    LC: QA297.T58.
c
c    Joan Westlake,
c    A Handbook of Numerical Matrix Inversion and Solution of 
c    Linear Equations,
c    John Wiley, 1968,
c    ISBN13: 978-0471936756,
c    LC: QA263.W47.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer MU, the number of superdiagonals.
c    MU must be at least 0, and no more than N-1.
c
c    Output, double precision A(MU+1,N), the matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      integer i
      integer j
      integer m

      do j = 1, n
        do i = 1, mu + 1
          a(i,j) = 0.0D+00
        end do
      end do

      do j = 2, n
        a(mu,j) = -1.0D+00
      end do

      do j = 1, n
        a(mu+1,j) =  +2.0D+00
      end do
   
      return
      end
      subroutine r8pbu_mv ( m, n, mu, a, x, b )

c*********************************************************************72
c
cc R8PBU_MV multiplies an R8PBU matrix by an R8VEC.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer MU, the number of superdiagonals in the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBU matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the result vector A * x.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(n)
      integer i
      integer ieqn
      integer j
      integer m
      double precision x(n)
c
c  Multiply X by the diagonal of the matrix.
c
      do j = 1, n
        b(j) = a(mu+1,j) * x(j)
      end do
c
c  Multiply X by the superdiagonals of the matrix.
c
      do i = mu, 1, -1
        do j = mu + 2 - i, n
          ieqn = i + j - mu - 1
          b(ieqn) = b(ieqn) + a(i,j) * x(j)
          b(j) = b(j) + a(i,j) * x(ieqn)
        end do
      end do

      return
      end
      subroutine r8pbu_resid ( m, n, mu, a, x, b, r )

c*********************************************************************72
c
cc R8PBU_RESID computes the residual R = B-A*X for R8PBU matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer MU, the number of superdiagonals in the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Input, double precision A(MU+1,N), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Input, double precision B(M), the desired result A * x.
c
c    Output, double precision R(M), the residual R = B - A * X.
c
      implicit none

      integer m
      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(m)
      integer i
      double precision r(m)
      double precision x(n)

      call r8pbu_mv ( m, n, mu, a, x, r )

      do i = 1, m
        r(i) = b(i) - r(i)
      end do

      return
      end
      subroutine r8sd_cg ( n, ndiag, offset, a, b, x )

c*********************************************************************72
c
cc R8SD_CG uses the conjugate gradient method on an R8SD linear system.
c
c  Discussion:
c
c    The R8SD storage format is for symmetric matrices whose only nonzero
c    entries occur along a few diagonals, but for which these diagonals are 
c    not all close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0, and 
c    each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c
c    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
c    we then create an array B that has N rows and NDIAG columns, and simply
c    "collapse" the matrix A to the left:
c
c    For the conjugate gradient method to be applicable, the matrix A must 
c    be a positive definite symmetric matrix.
c
c    The method is designed to reach the solution to the linear system
c      A * x = b
c    after N computational steps.  However, roundoff may introduce
c    unacceptably large errors for some problems.  In such a case,
c    calling the routine a second time, using the current solution estimate
c    as the new starting guess, should result in improved results.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Frank Beckman,
c    The Solution of Linear Equations by the Conjugate Gradient Method,
c    in Mathematical Methods for Digital Computers,
c    edited by John Ralston, Herbert Wilf,
c    Wiley, 1967,
c    ISBN: 0471706892,
c    LC: QA76.5.R3.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer NDIAG, the number of diagonals that are stored.
c    NDIAG must be at least 1 and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8SD matrix.
c
c    Input, double precision B(N), the right hand side vector.
c
c    Input/output, double precision X(N).
c    On input, an estimate for the solution, which may be 0.
c    On output, the approximate solution vector.  Note that repeated
c    calls to this routine, using the value of X output on the previous
c    call, MAY improve the solution.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      integer i
      integer it
      integer m
      integer offset(ndiag)
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      double precision r8vec_dot_product
      double precision rap
      double precision x(n)
c
c  Initialize
c    AP = A * x,
c    R  = b - A * x,
c    P  = b - A * x.
c
      call r8sd_mv ( n, n, ndiag, offset, a, x, ap )

      do i = 1, n 
        r(i) = b(i) - ap(i)
      end do

      do i = 1, n
        p(i) = b(i) - ap(i)
      end do
c
c  Do the N steps of the conjugate gradient method.
c
      do it = 1, n
c
c  Compute the matrix*vector product AP = A*P.
c
        call r8sd_mv ( n, n, ndiag, offset, a, p, ap )
c
c  Compute the dot products
c    PAP = P*AP,
c    PR  = P*R
c  Set
c    ALPHA = PR / PAP.
c
        pap = r8vec_dot_product ( n, p, ap )
        pr = r8vec_dot_product ( n, p, r )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
c
c  Set
c    X = X + ALPHA * P
c    R = R - ALPHA * AP.
c
        do i = 1, n
          x(i) = x(i) + alpha * p(i)
        end do

        do i = 1, n
          r(i) = r(i) - alpha * ap(i)
        end do
c
c  Compute the vector dot product
c    RAP = R*AP
c  Set
c    BETA = - RAP / PAP.
c
        rap = r8vec_dot_product ( n, r, ap )

        beta = - rap / pap
c
c  Update the perturbation vector
c    P = R + BETA * P.
c
        do i = 1, n
          p(i) = r(i) + beta * p(i)
        end do

      end do

      return
      end
      subroutine r8sd_dif2 ( m, n, ndiag, offset, a )

c*********************************************************************72
c
cc R8SD_DIF2 returns the DIF2 matrix in R8SD format.
c
c  Example:
c
c    N = 5
c
c    2 -1  .  .  .
c   -1  2 -1  .  .
c    . -1  2 -1  .
c    .  . -1  2 -1
c    .  .  . -1  2
c
c  Properties:
c
c    A is banded, with bandwidth 3.
c
c    A is tridiagonal.
c
c    Because A is tridiagonal, it has property A (bipartite).
c
c    A is a special case of the TRIS or tridiagonal scalar matrix.
c
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c
c    A is Toeplitz: constant along diagonals.
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    A is positive definite.
c
c    A is an M matrix.
c
c    A is weakly diagonally dominant, but not strictly diagonally dominant.
c
c    A has an LU factorization A = L * U, without pivoting.
c
c      The matrix L is lower bidiagonal with subdiagonal elements:
c
c        L(I+1,I) = -I/(I+1)
c
c      The matrix U is upper bidiagonal, with diagonal elements
c
c        U(I,I) = (I+1)/I
c
c      and superdiagonal elements which are all -1.
c
c    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
c
c      L(I,I) =    sqrt ( (I+1) / I )
c      L(I,I-1) = -sqrt ( (I-1) / I )
c
c    The eigenvalues are
c
c      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
c                = 4 SIN^2(I*PI/(2*N+2))
c
c    The corresponding eigenvector X(I) has entries
c
c       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
c
c    Simple linear systems:
c
c      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
c
c      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
c
c    det ( A ) = N + 1.
c
c    The value of the determinant can be seen by induction,
c    and expanding the determinant across the first row:
c
c      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
c                = 2 * N - (N-1)
c                = N + 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 July 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Gregory, David Karney,
c    A Collection of Matrices for Testing Computational Algorithms,
c    Wiley, 1969,
c    ISBN: 0882756494,
c    LC: QA263.68
c
c    Morris Newman, John Todd,
c    Example A8,
c    The evaluation of matrix inversion programs,
c    Journal of the Society for Industrial and Applied Mathematics,
c    Volume 6, Number 4, pages 466-476, 1958.
c
c    John Todd,
c    Basic Numerical Mathematics,
c    Volume 2: Numerical Algebra,
c    Birkhauser, 1980,
c    ISBN: 0817608117,
c    LC: QA297.T58.
c
c    Joan Westlake,
c    A Handbook of Numerical Matrix Inversion and Solution of 
c    Linear Equations,
c    John Wiley, 1968,
c    ISBN13: 978-0471936756,
c    LC: QA263.W47.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer NDIAG, the number of diagonals that are stored.
c    NDIAG must be at least 1 and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.  It is simply assumed that OFFSET(1) = 0 and OFFSET(2) = 1.
c
c    Output, double precision A(N,NDIAG), the R8SD matrix.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      integer i
      integer j
      integer m
      integer offset(ndiag)

      do j = 1, ndiag
        do i = 1, n
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        a(i,1) =  2.0D+00
      end do

      do i = 1, n - 1
        a(i,2) = -1.0D+00
      end do
   
      return
      end
      subroutine r8sd_mv ( m, n, ndiag, offset, a, x, b )

c*********************************************************************72
c
cc R8SD_MV multiplies an R8SD matrix by an R8VEC.
c
c  Discussion:
c
c    The R8SD storage format is for symmetric matrices whose only nonzero 
c    entries occur along a few diagonals, but for which these diagonals are not 
c    all close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0, and 
c    each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c
c    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
c    we then create an array B that has N rows and NDIAG columns, and simply
c    "collapse" the matrix A to the left:
c
c  Example:
c
c    The "offset" value is printed above each column.
c
c    Original matrix               New Matrix
c
c       0   1   2   3   4   5       0   1   3   5
c
c      11  12   0  14   0  16      11  12  14  16
c      21  22  23   0  25   0      22  23  25  --
c       0  32  33  34   0  36      33  34  36  --
c      41   0  43  44  45   0      44  45  --  --
c       0  52   0  54  55  56      55  56  --  --
c      61   0  63   0  65  66      66  --  --  --
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer NDIAG, the number of diagonals that are stored.
c    NDIAG must be at least 1 and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8SD matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer m
      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision b(n)
      integer i
      integer j
      integer jdiag
      integer offset(ndiag)
      double precision x(n)

      do i = 1, n
        b(i) = 0.0D+00
      end do

      do i = 1, n
        do jdiag = 1, ndiag
          if ( 0 .le. offset(jdiag) ) then
            j = i + offset(jdiag)
            if ( 1 .le. j .and. j .le. n ) then
              b(i) = b(i) + a(i,jdiag) * x(j)
              if ( offset(jdiag) .ne. 0 ) then
                b(j) = b(j) + a(i,jdiag) * x(i)
              end if
            end if
          end if
        end do
      end do

      return
      end
      subroutine r8sd_resid ( m, n, ndiag, offset, a, x, b, r )

c*********************************************************************72
c
cc R8SD_RESID computes the residual R = B-A*X for R8SD matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer NDIAG, the number of diagonals that are stored.
c    NDIAG must be at least 1 and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8SD matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Input, double precision B(M), the desired result A * x.
c
c    Output, double precision R(M), the residual R = B - A * X.
c
      implicit none

      integer m
      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision b(m)
      integer i
      integer offset(ndiag)
      double precision r(m)
      double precision x(n)

      call r8sd_mv ( m, n, ndiag, offset, a, x, r )

      do i = 1, m
        r(i) = b(i) - r(i)
      end do

      return
      end
      subroutine r8sp_cg ( n, nz_num, row, col, a, b, x )

c*********************************************************************72
c
cc R8SP_CG uses the conjugate gradient method on an R8SP system.
c
c  Discussion:
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c    The matrix A must be a positive definite symmetric band matrix.
c
c    The method is designed to reach the solution after N computational
c    steps.  However, roundoff may introduce unacceptably large errors for
c    some problems.  In such a case, calling the routine again, using
c    the computed solution as the new starting estimate, should improve
c    the results.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Frank Beckman,
c    The Solution of Linear Equations by the Conjugate Gradient Method,
c    in Mathematical Methods for Digital Computers,
c    edited by John Ralston, Herbert Wilf,
c    Wiley, 1967,
c    ISBN: 0471706892,
c    LC: QA76.5.R3.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer NZ_NUM, the number of nonzero elements in
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, double precision B(N), the right hand side vector.
c
c    Input/output, double precision X(N).
c    On input, an estimate for the solution, which may be 0.
c    On output, the approximate solution vector.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      integer col(nz_num)
      integer i
      integer it
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      double precision r8vec_dot_product
      integer row(nz_num)
      double precision rap
      double precision x(n)
c
c  Initialize
c    AP = A * x,
c    R  = b - A * x,
c    P  = b - A * x.
c
      call r8sp_mv ( n, n, nz_num, row, col, a, x, ap )

      do i = 1, n 
        r(i) = b(i) - ap(i)
      end do

      do i = 1, n
        p(i) = b(i) - ap(i)
      end do
c
c  Do the N steps of the conjugate gradient method.
c
      do it = 1, n
c
c  Compute the matrix*vector product AP=A*P.
c
        call r8sp_mv ( n, n, nz_num, row, col, a, p, ap )
c
c  Compute the dot products
c    PAP = P*AP,
c    PR  = P*R
c  Set
c    ALPHA = PR / PAP.
c
        pap = r8vec_dot_product ( n, p, ap )
        pr = r8vec_dot_product ( n, p, r )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
c
c  Set
c    X = X + ALPHA * P
c    R = R - ALPHA * AP.
c
        do i = 1, n
          x(i) = x(i) + alpha * p(i)
        end do

        do i = 1, n
          r(i) = r(i) - alpha * ap(i)
        end do
c
c  Compute the vector dot product
c    RAP = R*AP
c  Set
c    BETA = - RAP / PAP.
c
        rap = r8vec_dot_product ( n, r, ap )

        beta = - rap / pap
c
c  Update the perturbation vector
c    P = R + BETA * P.
c
        do i = 1, n
          p(i) = r(i) + beta * p(i)
        end do

      end do

      return
      end
      subroutine r8sp_dif2 ( m, n, nz_num, row, col, a )

c*********************************************************************72
c
cc R8SP_DIF2 returns the DIF2 matrix in R8SP format.
c
c  Example:
c
c    N = 5
c
c    2 -1  .  .  .
c   -1  2 -1  .  .
c    . -1  2 -1  .
c    .  . -1  2 -1
c    .  .  . -1  2
c
c  Properties:
c
c    A is banded, with bandwidth 3.
c
c    A is tridiagonal.
c
c    Because A is tridiagonal, it has property A (bipartite).
c
c    A is a special case of the TRIS or tridiagonal scalar matrix.
c
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c
c    A is Toeplitz: constant along diagonals.
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    A is positive definite.
c
c    A is an M matrix.
c
c    A is weakly diagonally dominant, but not strictly diagonally dominant.
c
c    A has an LU factorization A = L * U, without pivoting.
c
c      The matrix L is lower bidiagonal with subdiagonal elements:
c
c        L(I+1,I) = -I/(I+1)
c
c      The matrix U is upper bidiagonal, with diagonal elements
c
c        U(I,I) = (I+1)/I
c
c      and superdiagonal elements which are all -1.
c
c    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
c
c      L(I,I) =    sqrt ( (I+1) / I )
c      L(I,I-1) = -sqrt ( (I-1) / I )
c
c    The eigenvalues are
c
c      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
c                = 4 SIN^2(I*PI/(2*N+2))
c
c    The corresponding eigenvector X(I) has entries
c
c       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
c
c    Simple linear systems:
c
c      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
c
c      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
c
c    det ( A ) = N + 1.
c
c    The value of the determinant can be seen by induction,
c    and expanding the determinant across the first row:
c
c      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
c                = 2 * N - (N-1)
c                = N + 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 July 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Gregory, David Karney,
c    A Collection of Matrices for Testing Computational Algorithms,
c    Wiley, 1969,
c    ISBN: 0882756494,
c    LC: QA263.68
c
c    Morris Newman, John Todd,
c    Example A8,
c    The evaluation of matrix inversion programs,
c    Journal of the Society for Industrial and Applied Mathematics,
c    Volume 6, Number 4, pages 466-476, 1958.
c
c    John Todd,
c    Basic Numerical Mathematics,
c    Volume 2: Numerical Algebra,
c    Birkhauser, 1980,
c    ISBN: 0817608117,
c    LC: QA297.T58.
c
c    Joan Westlake,
c    A Handbook of Numerical Matrix Inversion and Solution of 
c    Linear Equations,
c    John Wiley, 1968,
c    ISBN13: 978-0471936756,
c    LC: QA263.W47.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer NZ_NUM, the number of nonzero elements in
c    the matrix.
c
c    Output, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Output, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer i
      integer k
      integer m
      integer row(nz_num)

      k = 0
      do i = 1, m

        if ( 0 .lt. i - 1 ) then
          k = k + 1
          row(k) = i
          col(k) = i - 1
          a(k) = -1.0D+00
        end if

        k = k + 1
        row(k) = i
        col(k) = i
        a(k) = 2.0D+00

        if ( i .lt. n ) then
          k = k + 1
          row(k) = i
          col(k) = i + 1
          a(k) = -1.0D+00
        end if

      end do

      return
      end
      subroutine r8sp_mv ( m, n, nz_num, row, col, a, x, b )

c*********************************************************************72
c
cc R8SP_MV multiplies an R8SP matrix by an R8VEC.
c
c  Discussion:
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product vector A*X.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(m)
      integer col(nz_num)
      integer i
      integer j
      integer k
      integer row(nz_num)
      double precision x(n)

      do i = 1, m
        b(i) = 0.0D+00
      end do

      do k = 1, nz_num

        i = row(k)
        j = col(k)
        b(i) = b(i) + a(k) * x(j)

      end do

      return
      end
      subroutine r8sp_resid ( m, n, nz_num, row, col, a, x, b, r )

c*********************************************************************72
c
cc R8SP_RESID computes the residual R = B-A*X for R8SP matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2014
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer NZ_NUM, the number of nonzero elements in
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Input, double precision B(M), the desired result A * x.
c
c    Output, double precision R(M), the residual R = B - A * X.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(m)
      integer col(nz_num)
      integer i
      double precision r(m)
      integer row(nz_num)
      double precision x(n)

      call r8sp_mv ( m, n, nz_num, row, col, a, x, r )

      do i = 1, m
        r(i) = b(i) - r(i)
      end do

      return
      end
      function r8vec_diff_norm ( n, a, b )

c*********************************************************************72
c
cc R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    The vector L2 norm is defined as:
c
c      R8VEC_NORM_L2 = sqrt ( sum ( 1 .le. I .le. N ) A(I)^2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, double precision A(N), B(N), the vectors.
c
c    Output, double precision R8VEC_DIFF_NORM, the L2 norm of A - B.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      double precision r8vec_diff_norm
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + ( a(i) - b(i) )**2
      end do
      value = sqrt ( value )

      r8vec_diff_norm = value

      return
      end
      function r8vec_dot_product ( n, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    In FORTRAN90, the system routine DOT_PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), V2(N), the vectors.
c
c    Output, double precision R8VEC_DOT_PRODUCT, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot_product
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r8vec_dot_product = value

      return
      end
      subroutine r8vec_house_column ( n, a, k, v )

c*********************************************************************72
c
cc R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
c
c  Discussion:
c
c    The routine returns a vector V that defines a Householder
c    premultiplier matrix H(V) that zeros out the subdiagonal entries of
c    column K of the matrix A.
c
c       H(V) = I - 2 * v * v'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix A.
c
c    Input, double precision A(N), column K of the matrix A.
c
c    Input, integer K, the column of the matrix to be modified.
c
c    Output, double precision V(N), a vector of unit L2 norm which defines an
c    orthogonal Householder premultiplier matrix H with the property
c    that the K-th column of H*A is zero below the diagonal.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer k
      double precision s
      double precision v(n)
      double precision vnorm

      do i = 1, n
        v(i) = 0.0D+00
      end do

      if ( k .lt. 1 .or. n .le. k ) then
        return
      end if

      s = 0.0D+00
      do i = k, n
        s = s + a(i)**2
      end do
      s = sqrt ( s )

      if ( s .eq. 0.0D+00 ) then
        return
      end if

      v(k) = a(k) + sign ( s, a(k) )
      do i = k + 1, n
        v(i) = a(i)
      end do

      vnorm = 0.0D+00
      do i = k, n
        vnorm = vnorm + v(i) * v(i)
      end do
      vnorm = sqrt ( vnorm )

      do i = k, n
        v(i) = v(i) / vnorm
      end do

      return
      end
      function r8vec_norm ( n, a )

c*********************************************************************72
c
cc R8VEC_NORM returns the L2 norm of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    The vector L2 norm is defined as:
c
c      R8VEC_NORM = sqrt ( sum ( 1 .le. I .le. N ) A(I)^2 ).
c
c  Licensing:
c
c    This code is distributed under the MIT license.
c
c  Modified:
c
c    BY Sourangshu Ghosh
c
c  Author:
c
c  Sourangshu Ghosh
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, double precision A(N), the vector whose L2 norm is desired.
c
c    Output, double precision R8VEC_NORM, the L2 norm of A.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision r8vec_norm
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + a(i) * a(i)
      end do
      value = sqrt ( value )

      r8vec_norm = value

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the MIT license.
c
c  Modified:
c
c    By Sourangshu Ghosh
c
c  Author:
c
c    Sourangshu Ghosh
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the MIT license.
c
c  Modified:
c
c    By Sourangshu Ghosh
c
c  Author:
c
c    Sourangshu Ghosh
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer k
      integer seed
      double precision r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the MIT license.
c
c  Modified:
c
c    By Sourangshu Ghosh
c
c  Author:
c
c    Sourangshu Ghosh
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
