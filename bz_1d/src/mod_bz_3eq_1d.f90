module mod_bz_3eq_1d

use mod_precision
use mod_cartesian_grid

integer, private, parameter :: neq = 3
real(kind=dp), private, parameter :: mu = 1.d-05
real(kind=dp), private, parameter :: eps = 1.d-02
real(kind=dp), private, parameter :: f = 3.d0
real(kind=dp), private, parameter :: q = 2.d-04
real(kind=dp), private, parameter :: da = 1.d0/400.d0
real(kind=dp), private, parameter :: db = 1.d0/400.d0
real(kind=dp), private, parameter :: dc = 6.d-01/400.d0
real(kind=dp), private :: daoverdxdx
real(kind=dp), private :: dboverdxdx
real(kind=dp), private :: dcoverdxdx
real(kind=dp), private :: oneovermu
real(kind=dp), private :: oneovereps
integer, private :: nx
real(kind=dp), private :: xmin

contains

subroutine init_bz_3eq(grid)
  implicit none
  type(cartesian_grid_type) :: grid
  real(kind=dp) :: dx

  dx = grid%dx

  daoverdxdx = da / (dx*dx)
  dboverdxdx = db / (dx*dx)
  dcoverdxdx = dc / (dx*dx)

  oneovermu  = 1.d0 / mu
  oneovereps = 1.d0 / eps

  nx = grid%nx

  xmin = grid%xmin

end subroutine init_bz_3eq

function f_bz_3eq(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_bz_3eq
  ! local 
  real(kind=dp) :: aim1, ai, aip1
  real(kind=dp) :: bim1, bi, bip1
  real(kind=dp) :: cim1, ci, cip1
  integer :: inx, irow

  ! left boundary
  ai   = u(1)
  bi   = u(2)
  ci   = u(3)
  aip1 = u(4)
  bip1 = u(5)
  cip1 = u(6)

  f_bz_3eq(1) = (daoverdxdx * (aip1 - ai)) + oneovermu  * (-q*ai - ai*bi + f*ci)
  f_bz_3eq(2) = (dboverdxdx * (bip1 - bi)) + oneovereps * ( q*ai - ai*bi + bi - bi*bi)
  f_bz_3eq(3) = (dcoverdxdx * (cip1 - ci)) + (bi - ci)

  do inx = 2, nx-1

    irow = 1 + (inx-1)*neq

    aim1 = u(irow-3)
    bim1 = u(irow-2)
    cim1 = u(irow-1)
    ai   = u(irow)
    bi   = u(irow+1)
    ci   = u(irow+2)
    aip1 = u(irow+3)
    bip1 = u(irow+4)
    cip1 = u(irow+5)

    f_bz_3eq(irow)   = (daoverdxdx * (aim1 - 2.d0*ai + aip1)) + oneovermu  * (-q*ai - ai*bi + f*ci)
    f_bz_3eq(irow+1) = (dboverdxdx * (bim1 - 2.d0*bi + bip1)) + oneovereps * ( q*ai - ai*bi + bi - bi*bi)
    f_bz_3eq(irow+2) = (dcoverdxdx * (cim1 - 2.d0*ci + cip1)) + (bi - ci)

  end do 

  ! right boundary
  aim1 = u(n-5)
  bim1 = u(n-4)
  cim1 = u(n-3)
  ai   = u(n-2)
  bi   = u(n-1)
  ci   = u(n)

  f_bz_3eq(n-2) = (daoverdxdx * (aim1 - ai)) + oneovermu  * (-q*ai - ai*bi + f*ci)
  f_bz_3eq(n-1) = (dboverdxdx * (bim1 - bi)) + oneovereps * ( q*ai - ai*bi + bi - bi*bi)
  f_bz_3eq(n)   = (dcoverdxdx * (cim1 - ci)) + (bi - ci)

end function f_bz_3eq

function f_loc_bz_3eq_reac(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_loc_bz_3eq_reac
  ! local
  real(kind=dp) :: ai
  real(kind=dp) :: bi
  real(kind=dp) :: ci

  ai   = u(1)
  bi   = u(2)
  ci   = u(3)

  f_loc_bz_3eq_reac(1) = oneovermu  * (-q*ai - ai*bi + f*ci)
  f_loc_bz_3eq_reac(2) = oneovereps * ( q*ai - ai*bi + bi - bi*bi)
  f_loc_bz_3eq_reac(3) = (bi - ci)

end function f_loc_bz_3eq_reac

function f_bz_3eq_diff(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_bz_3eq_diff
  ! local 
  real(kind=dp) :: aim1, ai, aip1
  real(kind=dp) :: bim1, bi, bip1
  real(kind=dp) :: cim1, ci, cip1
  integer :: inx, irow

  ! left boundary
  ai   = u(1)
  bi   = u(2)
  ci   = u(3)
  aip1 = u(4)
  bip1 = u(5)
  cip1 = u(6)

  f_bz_3eq_diff(1) = (daoverdxdx * (aip1 - ai))
  f_bz_3eq_diff(2) = (dboverdxdx * (bip1 - bi))
  f_bz_3eq_diff(3) = (dcoverdxdx * (cip1 - ci))

  do inx = 2, nx-1

    irow = 1 + (inx-1)*neq

    aim1 = u(irow-3)
    bim1 = u(irow-2)
    cim1 = u(irow-1)
    ai   = u(irow)
    bi   = u(irow+1)
    ci   = u(irow+2)
    aip1 = u(irow+3)
    bip1 = u(irow+4)
    cip1 = u(irow+5)

    f_bz_3eq_diff(irow)   = (daoverdxdx * (aim1 - 2.d0*ai + aip1))
    f_bz_3eq_diff(irow+1) = (dboverdxdx * (bim1 - 2.d0*bi + bip1))
    f_bz_3eq_diff(irow+2) = (dcoverdxdx * (cim1 - 2.d0*ci + cip1))

  end do 

  ! right boundary
  aim1 = u(n-5)
  bim1 = u(n-4)
  cim1 = u(n-3)
  ai   = u(n-2)
  bi   = u(n-1)
  ci   = u(n)

  f_bz_3eq_diff(n-2) = (daoverdxdx * (aim1 - ai)) 
  f_bz_3eq_diff(n-1) = (dboverdxdx * (bim1 - bi))
  f_bz_3eq_diff(n)   = (dcoverdxdx * (cim1 - ci))

end function f_bz_3eq_diff

subroutine bz_3eq_init_sol(u)
  implicit none
  real(kind=dp), dimension(:) :: u
  ! local
  real(kind=dp), allocatable, dimension(:) :: a, b, c
  real(kind=dp) :: xcoor, ycoor, ylim
  real(kind=dp), parameter :: pi = 4.d0 * atan(1.d0)
  integer :: inx, irow

  allocate(a(nx))
  allocate(b(nx))
  allocate(c(nx))

  ylim = 0.05d0

  do  inx = 1, int(nx/20)
    xcoor = 0.5d0
    ycoor = dfloat(inx)/dfloat(nx/20) - ylim

    if (ycoor >= 0.d0 .and. ycoor<= 0.3*xcoor) then
      b(inx) = 0.8d0
    else
      b(inx) = q*(f+1.d0)/(f-1.d0)
    endif

   if ( ycoor>=0.d0 ) then
     c(inx) = q*(f+1.d0)/(f-1.d0) + datan(ycoor/xcoor)/(8.d0*pi*f)
   else
     c(inx) = q*(f+1.d0)/(f-1.d0) + (datan(ycoor/xcoor) + 2.d0*pi)/(8.d0*pi*f)
   endif

  enddo

  do inx=int(nx/20)+1, nx
    b(inx)=b(int(nx/20))
    c(inx)=c(int(nx/20))
  end do

  do inx=1, nx
    a(inx)=(f*c(inx))/(q+b(inx))
  end do

  do inx = 1, nx
    irow = 1 + (inx-1)*neq
    u(irow)   = a(inx)
    u(irow+1) = b(inx)
    u(irow+2) = c(inx)
  end do

  deallocate(a, b, c)

end subroutine bz_3eq_init_sol

end module mod_bz_3eq_1d
