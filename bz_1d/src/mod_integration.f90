module mod_integration

use mod_precision
use mod_radau
use mod_strang
use mod_strang_drd
use mod_lie_rd
use mod_lie_dr
use mod_bz_2eq_1d
use mod_bz_3eq_1d

contains

subroutine integrate(method, tol, neq, nx, u, nt, tini, tend)
  implicit none
  ! arguments
  character(len=20) :: method
  real(kind=dp) :: tol
  integer :: neq
  integer :: nx
  real(kind=dp), dimension(:) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  integer :: ntot

  ntot = neq * nx

  select case (method) 
    case ("radau5")
      print *
      print *, "Radau5 integration"
      if (neq == 2) then
        call radau5_integration(f_bz_2eq, tol, neq, ntot, u, tini, tend, .true.)
      else
        call radau5_integration(f_bz_3eq, tol, neq, ntot, u, tini, tend, .true.)
      end if
    case ("strang")
      print *
      print *, "Strang reaction-diffusion-diffusion splitting integration"
      print *, "  RADAU5 for reaction"
      print *, "  ROCK4  for diffusion"
      if (neq == 2) then
        call strang_integration(f_loc_bz_2eq_reac, f_bz_2eq_diff, tol, neq, nx, u, nt, tini, tend) 
      else
        call strang_integration(f_loc_bz_3eq_reac, f_bz_3eq_diff, tol, neq, nx, u, nt, tini, tend) 
      end if 
    case ("strang_drd")
      print *
      print *, "Strang diffusion-reaction-diffusion splitting integration"
      print *, "  RADAU5 for reaction"
      print *, "  ROCK4  for diffusion"
      if (neq == 2) then
        call strang_integration_drd(f_loc_bz_2eq_reac, f_bz_2eq_diff, tol, neq, nx, u, nt, tini, tend) 
      else
        call strang_integration_drd(f_loc_bz_3eq_reac, f_bz_3eq_diff, tol, neq, nx, u, nt, tini, tend) 
      end if 

    case ("lie_rd")
      print *
      print *, "Lie reaction-diffusion splitting integration"
      print *, "  RADAU5 for reaction"
      print *, "  ROCK4  for diffusion"
      if (neq == 2) then
        call lie_integration_rd(f_loc_bz_2eq_reac, f_bz_2eq_diff, tol, neq, nx, u, nt, tini, tend) 
      else
        call lie_integration_rd(f_loc_bz_3eq_reac, f_bz_3eq_diff, tol, neq, nx, u, nt, tini, tend) 
      end if
    case ("lie_dr")
      print *
      print *, "Lie diffusion-reaction splitting integration"
      print *, "  RADAU5 for reaction"
      print *, "  ROCK4  for diffusion"
      if (neq == 2) then
        call lie_integration_dr(f_loc_bz_2eq_reac, f_bz_2eq_diff, tol, neq, nx, u, nt, tini, tend) 
      else
        call lie_integration_dr(f_loc_bz_3eq_reac, f_bz_3eq_diff, tol, neq, nx, u, nt, tini, tend) 
      end if  
    case default
      print *
      print *, method, "Unknown integration method"
      call exit(0)
  end select

end subroutine integrate

end module mod_integration
