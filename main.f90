! 2D diffusion equation

module functions_module
  contains
    subroutine thomas_algorithm(a, b, c, r, x, s, e)
      real *8, allocatable, intent(in):: a(:), c(:)
      real *8, allocatable, intent(inout):: x(:), b(:), r(:)
      integer, intent(in):: s, e

      integer:: i

      do i=s+1, e
        b(i) = b(i) - a(i)/b(i-1)*c(i-1)
        r(i) = r(i) - a(i)/b(i-1)*r(i-1)
      end do

      x(e) = r(e)/b(e)
      do i=e-1,s,-1
        x(i) = (r(i)-c(i)*x(i+1))/b(i)
      end do

      return
    end subroutine thomas_algorithm


    subroutine save_numeric(nu, n)
      real *8, allocatable, intent(in):: nu(:,:)
      integer, intent(in):: n

      integer:: i, j

      open(unit=2, file="numeric.dat")

      do i=0, n
        do j=0, n
          write (2,*) i, j, nu(i, j)
        end do
      end do

      return
    end subroutine save_numeric
end module functions_module

module diffusion_module
  contains
    subroutine crank_nicolson(nu, n, s)
      use functions_module
      real *8, allocatable, intent(inout):: nu(:,:)
      real *8, intent(in):: s
      integer, intent(in):: n

      real *8, allocatable:: a(:), b(:), c(:), r(:), x(:)
      real *8, allocatable:: z(:,:)
      integer:: i, j

      allocate(z(0:n, 0:n))

      do j=0, n
        do i=0, n
          z(i, j) = 0.0d0
        end do
      end do

      do j=1, n-1
        allocate(a(1:n-1), b(1:n-1), c(1:n-1), r(1:n-1), x(1:n-1))

        do i=1, n-1
          a(i) = -s
          b(i) = (1.0d0 + 2.0d0*s)
          c(i) = -s
          r(i) = (s*nu(i, j-1)) + ((1.0d0 - 2.0d0*s)*nu(i, j)) + (s*nu(i, j+1))
        end do

        call thomas_algorithm(a, b, c, r, x, 1, n-1)

        z(0, j) = 0.0d0
        do i=1, n-1
          z(i, j) = x(i)
        end do
        z(n, j) = 0.0d0

        deallocate(a, b, c, r, x)
      end do

      do i=1, n-1
        allocate(a(1:n-1), b(1:n-1), c(1:n-1), r(1:n-1), x(1:n-1))

        do j=1, n-1
          a(j) = -s
          b(j) = (1.0d0 + 2.0d0*s)
          c(j) = -s
          r(j) = (s*z(i-1, j)) + ((1.0d0 - 2.0d0*s)*z(i, j)) + (s*z(i+1, j))
        end do

        call thomas_algorithm(a, b, c, r, x, 1, n-1)

        nu(i, 0) = 0.0d0
        do j=1, n-1
          nu(i, j) = x(j)
        end do
        nu(i, n) = 0.0d0

        deallocate(a, b, c, r, x)
      end do

      deallocate(z)

      return
    end subroutine crank_nicolson
end module diffusion_module

program diffusion_advection
  use functions_module
  use diffusion_module

  integer:: n, o, k, st, i, j
  real *8:: d, dt, l, s
  real *8, allocatable:: nu(:,:)

  l = 1.0d0
  st = 100

  do o=2, 6
    print *, "Execucao com o=", o
    n = 2 ** o
    d = l / (dfloat(n))
    dt = 0.0001d0
    s = dt / (2.0d0 * (d * d))

    allocate(nu(0:n, 0:n))

    !call initial_value(nu, n, d)

    do k=1, st-1
      call save_numeric(nu, n)

      call crank_nicolson(nu, n, s)
    end do

    deallocate(nu)
  end do

end program diffusion_advection


!###############################################################################
real *8 function analytic_solution(x, y, t)
  real *8, intent(in):: x, y, t
  real *8:: pi = ACOS(DBLE(-1.0))

  analytic_solution = SIN(pi*x) * SIN(pi*y) * EXP(-2.0d0*(pi*pi)*t)
end function analytic_solution
