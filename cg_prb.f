      program main

c*********************************************************************72
c
cc MAIN is the main program for CG_PRB.
c
c  Discussion:
c
c    CG_PRB tests the CG library.
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
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CG_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the CG library.'

      call test01 ( )
      call test02 ( )
      call test023 ( )
      call test025 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CG_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests R8GE_CG for a full storage matrix.
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
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(n,n)
      double precision b(n)
      double precision e_norm
      integer i
      double precision r(n)
      double precision r_norm
      double precision r8vec_diff_norm
      double precision r8vec_norm
      integer seed
      double precision x1(n)
      double precision x2(n)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  Test R8GE_CG, applying CG to a full storage matrix.'
c
c  Choose a random positive definite symmetric matrix A.
c
      seed = 123456789
      call pds_random ( n, seed, a )
c
c  Choose a random solution.
c
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x1 )
c
c  Compute the corresponding right hand side.
c
      call r8ge_mv ( n, n, a, x1, b )
c
c  Call the CG routine.
c
      do i = 1, n
        x2(i) = 1.0D+00
      end do

      call r8ge_cg ( n, a, b, x2 )
c
c  Compute the residual.
c
      call r8ge_resid ( n, n, a, x2, b, r )
      r_norm = r8vec_norm ( n, r )
c
c  Compute the error.
c
      e_norm = r8vec_diff_norm ( n, x1, x2 )
c
c  Report.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of variables N = ', n
      write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
      write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R83_CG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2014
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(3,n)
      double precision b(n)
      double precision e_norm
      integer i
      double precision r(n)
      double precision r_norm
      double precision r8vec_diff_norm
      double precision r8vec_norm
      integer seed
      double precision x1(n)
      double precision x2(n)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test R83_CG, applying CG to an R83 matrix.'
c
c  Let A be the -1 2 -1 matrix.
c
      call r83_dif2 ( n, n, a )
c
c  Choose a random solution.
c
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x1 )
c
c  Compute the corresponding right hand side.
c
      call r83_mv ( n, n, a, x1, b )
c
c  Call the CG routine.
c
      do i = 1, n
        x2(i) = 1.0D+00
      end do

      call r83_cg ( n, a, b, x2 )
c
c  Compute the residual.
c
      call r83_resid ( n, n, a, x2, b, r )
      r_norm = r8vec_norm ( n, r )
c
c  Compute the error.
c
      e_norm = r8vec_diff_norm ( n, x1, x2 )
c
c  Report.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of variables N = ', n
      write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
      write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm

      return
      end
      subroutine test023 ( )

c*********************************************************************72
c
cc TEST023 tests R83S_CG.
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
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(3)
      double precision b(n)
      double precision e_norm
      integer i
      double precision r(n)
      double precision r_norm
      double precision r8vec_diff_norm
      double precision r8vec_norm
      integer seed
      double precision x1(n)
      double precision x2(n)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST023'
      write ( *, '(a)' ) 
     &  '  Test R83S_CG, applying CG to an R83S matrix.'
c
c  Let A be the -1 2 -1 matrix.
c
      call r83s_dif2 ( n, n, a )
c
c  Choose a random solution.
c
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x1 )
c
c  Compute the corresponding right hand side.
c
      call r83s_mv ( n, n, a, x1, b )
c
c  Call the CG routine.
c
      do i = 1, n
        x2(i) = 1.0D+00
      end do

      call r83s_cg ( n, a, b, x2 )
c
c  Compute the residual.
c
      call r83s_resid ( n, n, a, x2, b, r )
      r_norm = r8vec_norm ( n, r )
c
c  Compute the error.
c
      e_norm = r8vec_diff_norm ( n, x1, x2 )
c
c  Report.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of variables N = ', n
      write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
      write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm

      return
      end
      subroutine test025 ( )

c*********************************************************************72
c
cc TEST025 tests R83T_CG.
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
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(3,n)
      double precision b(n)
      double precision e_norm
      integer i
      double precision r(n)
      double precision r_norm
      double precision r8vec_diff_norm
      double precision r8vec_norm
      integer seed
      double precision x1(n)
      double precision x2(n)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST025'
      write ( *, '(a)' ) 
     &  '  Test R83T_CG, applying CG to an R83T matrix.'
c
c  Let A be the -1 2 -1 matrix.
c
      call r83t_dif2 ( n, n, a )
c
c  Choose a random solution.
c
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x1 )
c
c  Compute the corresponding right hand side.
c
      call r83t_mv ( n, n, a, x1, b )
c
c  Call the CG routine.
c
      do i = 1, n
        x2(i) = 1.0D+00
      end do

      call r83t_cg ( n, a, b, x2 )
c
c  Compute the residual.
c
      call r83t_resid ( n, n, a, x2, b, r )
      r_norm = r8vec_norm ( n, r )
c
c  Compute the error.
c
      e_norm = r8vec_diff_norm ( n, x1, x2 )
c
c  Report.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of variables N = ', n
      write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
      write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests R8PBU_CG.
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
      implicit none

      integer mu
      parameter ( mu = 1 )
      integer n
      parameter ( n = 10 )

      double precision a(mu+1,n)
      double precision b(n)
      double precision e_norm
      integer i
      double precision r(n)
      double precision r_norm
      double precision r8vec_diff_norm
      double precision r8vec_norm
      integer seed
      double precision x1(n)
      double precision x2(n)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  Test R8PBU_CG, applying CG to an R8PBU matrix.'
c
c  Let A be the -1 2 -1 matrix.
c
      call r8pbu_dif2 ( n, n, mu, a )
c
c  Choose a random solution.
c
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x1 )
c
c  Compute the corresponding right hand side.
c
      call r8pbu_mv ( n, n, mu, a, x1, b )
c
c  Call the CG routine.
c
      do i = 1, n
        x2(i) = 1.0D+00
      end do

      call r8pbu_cg ( n, mu, a, b, x2 )
c
c  Compute the residual.
c
      call r8pbu_resid ( n, n, mu, a, x2, b, r )
      r_norm = r8vec_norm ( n, r )
c
c  Compute the error.
c
      e_norm = r8vec_diff_norm ( n, x1, x2 )
c
c  Report.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of variables N = ', n
      write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
      write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests R8SD_CG.
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
      implicit none

      integer n
      parameter ( n = 10 )
      integer ndiag
      parameter ( ndiag = 2 )

      double precision a(n,ndiag)
      double precision b(n)
      double precision e_norm
      integer i
      integer offset(ndiag)
      double precision r(n)
      double precision r_norm
      double precision r8vec_diff_norm
      double precision r8vec_norm
      integer seed
      double precision x1(n)
      double precision x2(n)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  Test R8SD_CG, applying CG to an R8SD matrix.'

      offset(1) = 0
      offset(2) = 1
c
c  Set A to the [-1 2 -1] matrix.
c
      call r8sd_dif2 ( n, n, ndiag, offset, a )
c
c  Choose a random solution.
c
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x1 )
c
c  Compute the corresponding right hand side.
c
      call r8sd_mv ( n, n, ndiag, offset, a, x1, b )
c
c  Call the CG routine.
c
      do i = 1, n
        x2(i) = 1.0D+00
      end do

      call r8sd_cg ( n, ndiag, offset, a, b, x2 )
c
c  Compute the residual.
c
      call r8sd_resid ( n, n, ndiag, offset, a, x2, b, r )
      r_norm = r8vec_norm ( n, r )
c
c  Compute the error.
c
      e_norm = r8vec_diff_norm ( n, x1, x2 )
c
c  Report.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of variables N = ', n
      write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
      write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests R8SP_CG.
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
      implicit none

      integer n
      parameter ( n = 10 )
      integer nz_num
      parameter ( nz_num = 3 * n - 2 )

      double precision a(nz_num)
      double precision b(n)
      integer col(nz_num)
      double precision e_norm
      integer i
      double precision r(n)
      integer row(nz_num)
      double precision r_norm
      double precision r8vec_diff_norm
      double precision r8vec_norm
      integer seed
      double precision x1(n)
      double precision x2(n)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) 
     &  '  Test R8SP_CG, applying CG to an R8SP matrix.'
c
c  Set A to the [-1 2 -1] matrix.
c
      call r8sp_dif2 ( n, n, nz_num, row, col, a )
c
c  Choose a random solution.
c
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x1 )
c
c  Compute the corresponding right hand side.
c
      call r8sp_mv ( n, n, nz_num, row, col, a, x1, b )
c
c  Call the CG routine.
c
      do i = 1, n
        x2(i) = 1.0D+00
      end do

      call r8sp_cg ( n, nz_num, row, col, a, b, x2 )
c
c  Compute the residual.
c
      call r8sp_resid ( n, n, nz_num, row, col, a, x2, b, r )
      r_norm = r8vec_norm ( n, r )
c
c  Compute the error.
c
      e_norm = r8vec_diff_norm ( n, x1, x2 )
c
c  Report.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of variables N = ', n
      write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
      write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm

      return
      end
