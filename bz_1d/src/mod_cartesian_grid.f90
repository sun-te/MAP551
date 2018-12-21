module mod_cartesian_grid

type cartesian_grid_type
  ! minimal and maximal coordinates of domain
  real(kind=kind(0.d0)) :: xmin, xmax
  ! number of interior points 
  integer :: nx
  ! step size 
  real(kind=kind(0.d0)) :: dx
end type cartesian_grid_type

contains

subroutine init_cartesian_grid(xmin, xmax, nxib, grid)
  implicit none
  real(kind=kind(0.d0)) :: xmin, xmax
  integer :: nxib
  type(cartesian_grid_type) :: grid

  grid%xmin = xmin
  grid%xmax = xmax

  grid%nx = nxib-2

  grid%dx = (xmax-xmin) / (nxib-1)

end subroutine init_cartesian_grid

end module mod_cartesian_grid
