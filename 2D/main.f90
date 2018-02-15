! 2D diffusion-advection equation

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


    subroutine initial_value(f, n, h)
      real *8, allocatable, intent(inout):: f(:,:)
      real *8, intent(in):: h
      integer, intent(in):: n

      integer:: i, j

      do i=1, n
        do j=1, n
          f(i, j) = analytic_solution(((j-1)*h), ((i-1)*h), 0.0d0)
        end do
      end do

      return
    end subroutine initial_value

    subroutine save_numeric(nu, n, h)
      real *8, allocatable, intent(in):: nu(:,:)
      real *8, intent(in):: h
      integer, intent(in):: n

      integer:: i, j
      
      do i=1, n
        do j=1, n
          write (2,*) (i-1)*h, (j-1)*h, nu(i, j)
        end do
      end do
      
      return
    end subroutine save_numeric


    subroutine crank_nicolson(nu, n, s, p, q)
      real *8, allocatable, intent(inout):: nu(:,:)
      real *8, intent(in):: s, p, q
      integer, intent(in):: n

      real *8, allocatable:: a(:), b(:), c(:), r(:), x(:)
      real *8, allocatable:: z(:,:)
      real *8:: ua, ub
      integer:: i, j

      allocate(z(1:n, 1:n))

      do j=1, n
        do i=1, n
          z(i, j) = 0.0d0
        end do
      end do

      do j=2, n-1
        allocate(a(1:n), b(1:n), c(1:n), r(1:n), x(1:n))

        do i=1, n-1
          a(i) = -(p + .5d0*s)
          b(i) = (1.0d0 + s)
          c(i) = (p - .5d0*s)
        end do

        ua = nu(1, j) - (p+.5d0*s)*nu(1,j-1) + s*nu(1, j) - (.5d0*s-p)*nu(1, j+1)
        r(1)   = r(1) - a(1)*ua
        do i=2, n-1
          r(i) = ((p+.5d0*s)*nu(i, j-1)) + ((1.0d0 - s)*nu(i, j)) + ((.5d0*s-p)*nu(i, j+1))
        end do
        ub = nu(n, j) - (p+.5d0*s)*nu(n,j-1) + s*nu(n, j) - (.5d0*s-p)*nu(n, j+1)
        r(n) = r(n) - c(n)*ub 
        
        call thomas_algorithm(a, b, c, r, x, 2, n-1)

!        z(1, j) = 0.0d0
        do i=1, n
          z(i, j) = x(i)
        end do
!        z(n, j) = 0.0d0

        deallocate(a, b, c, r, x)
      end do

      do i=2, n-1
        allocate(a(1:n), b(1:n), c(1:n), r(1:n), x(1:n))

        do j=1, n-1
          a(j) = -(q + .5d0*s)
          b(j) = (1.0d0 + s)
          c(j) = (q - .5d0*s)
        end do

        ua = nu(i, 1) - (q+.5d0*s)*nu(i-1, 1) + s*nu(i, 1) - (.5d0*s-q)*nu(i+1, 1)
        r(1) = r(1) - a(1)*ua
        do j=2, n-1
          r(j) = ((q+.5d0*s)*z(i-1, j)) + ((1.0d0 - s)*z(i, j)) + ((.5d0*s-q)*z(i+1, j))
        end do
        ub = nu(i, n) - (q+.5d0*s)*nu(i-1, n) + s*nu(i, n) - (.5d0*s-q)*nu(i+1, n)
        r(n-1) = r(n-1) + c(n-1)*ub

        call thomas_algorithm(a, b, c, r, x, 2, n-1)

!        nu(i, 1) = 0.0d0
        do j=1, n
          nu(i, j) = x(j)
        end do
!        nu(i, n) = 0.0d0

        deallocate(a, b, c, r, x)
      end do

      deallocate(z)

      return
    end subroutine crank_nicolson

end module functions_module

program diffusion
  use functions_module

  integer:: n, o, k, st, i, j
  real *8:: d, dt, l, s, u, v, p, q
  real *8, allocatable:: nu(:,:)

  l = 1.0d0
  st = 125
  
  open(unit=2, file="numeric.dat")
  n = 2 ** 5
  d = l / (n-1)
  pe = 7.0d0
  u = 0.0d0
  v = 0.0d0
  dt = 0.001d0
  s = dt / (d * d)
  p = (pe*u*dt)/(4.0d0*d)
  q = (pe*v*dt)/(4.0d0*d)

  allocate(nu(1:n, 1:n))
  call initial_value(nu, n, d)

  do k=1, st-1
    call save_numeric(nu, n, d)

    call crank_nicolson(nu, n, s, p, q)
  end do

  deallocate(nu)
  close(2)

end program diffusion


!###############################################################################
real *8 function analytic_solution(x, y, t)
  real *8, intent(in):: x, y, t
  real *8:: pi = ACOS(DBLE(-1.0))

  analytic_solution = SIN(pi*x) * SIN(pi*y) * EXP(-2.0d0*(pi*pi)*t)
end function analytic_solution

