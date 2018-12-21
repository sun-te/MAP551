module mod_strang

use mod_precision
use mod_radau
use mod_rock

contains


subroutine strang_integration(f_loc_reac, f_diff, tol, neq, nx, u, nt, tini, tend)
  implicit none
  ! arguments
  interface
    function f_loc_reac(n,y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f_loc_reac
    end function f_loc_reac
    function f_diff(n, y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f_diff
    end function f_diff
  end interface
  real(kind=dp) :: tol
  integer :: neq
  integer :: nx
  real(kind=dp), dimension(:) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  ! local variables
  real(kind=dp) :: dt, t
  real(kind=dp), dimension(neq) :: uloc
  integer :: ntot
  integer :: it, inx, irow

  ntot = neq*nx

  dt = (tend-tini) / (nt-1)
  print '("   Splitting time step = ", e10.4)', dt

  do it = 1, nt-1

    t = tini + (it-1)*dt
    print '("   Iteration ", i5, "  t = ", e10.4)', it, t+dt
 
    !!print *, " Radau5 integration for reaction on dt/2"
    do inx = 1, nx

      irow = 1 + (inx-1)*neq
      
      !!uloc(1) = u(irow)  
      !!uloc(2) = u(irow+1)  
      uloc(:) = u(irow:irow+neq-1)

      call radau5_integration(f_loc_reac, tol, neq, neq, uloc, t, t+(dt/2), .false.)

      !!u(irow)   = uloc(1)
      !!u(irow+1) = uloc(2)  
      u(irow:irow+neq-1) = uloc(:)

    end do

    !!print *, " Rock4 integration for diffusion on dt"
    call rock4_integration(f_diff, tol, ntot, u, 1, t, t+dt, .false.)

    !!print *, " Radau5 integration for reaction on dt/2"
    do inx = 1, nx

      irow = 1 + (inx-1)*neq
      
      !!uloc(1) = u(irow)  
      !!uloc(2) = u(irow+1)  
      uloc(:) = u(irow:irow+neq-1)

      call radau5_integration(f_loc_reac, tol, neq, neq, uloc, t+(dt/2), t+dt, .false.)

      !!u(irow)   = uloc(1)
      !!u(irow+1) = uloc(2)  
      u(irow:irow+neq-1) = uloc(:)

    end do

  end do

end subroutine strang_integration

end module mod_strang
