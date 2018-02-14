module routines
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

    subroutine initial_condition(nu, n, dx)
      real *8, allocatable, intent(inout):: nu(:)
      real *8, intent(in):: dx
      integer, intent(in):: n

      integer:: i
      real *8:: pi

      pi = 4.0d0*DATAN(1.0d0)

      do i=1, n-1
        nu(i) = SIN(pi*i*dx)
      end do

      return
    end subroutine initial_condition

    subroutine save_numeric(nu, n)
      real *8, allocatable, intent(in):: nu(:)
      integer, intent(in):: n

      integer:: i
      do i=1, n-1
        write (2,*) i, nu(i)
      end do

      return
    end subroutine save_numeric

end module routines


program adv_dif
    use routines

    real *8, allocatable:: a(:), b(:), c(:), nu(:), x(:), y(:)
    real *8:: u, dx, dt, r, s, l
    integer:: i, o, ts, n


    ts = 150
    l = 1.0d0
    u = -1.0d0
    pe = 5.0d0
    n = 2 ** 8
    dx = l / (n-1)
    dt = 0.001d0
    r = (pe*u*dt) / (4.0d0*dx)
    s = dt / (dx*dx)

    open(unit=2, file='numeric.dat')
    allocate(nu(1:n-1))
    call initial_condition(nu, n, dx)

    do o=2, ts+1
      allocate(a(1:n-1), b(1:n-1), c(1:n-1), x(1:n-1), y(1:n-1))
      call save_numeric(nu, n)

      y(1) = 0.0d0
      do i=2, n-2
        a = -(r + .5d0*s)
        b = 1.0d0 + s
        c = r - .5d0*s
        y(i) = (r+.5d0*s)*nu(i-1) + (1.0d0 - s)*nu(i) + (.5d0*s - r)*nu(i+1)
      end do
      y(n-1) = 0.0d0

      call thomas_algorithm(a, b, c, y, x, 1, n-1)

      do i=1, n-1
        nu(i) = x(i)
      end do

      deallocate(a, b, c, x, y)
    end do
    deallocate(nu)
    close(2)

end program adv_dif

