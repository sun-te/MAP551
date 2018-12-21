module mod_radau

use mod_precision

abstract interface
  function func(n, y)
     integer :: n
     real(kind=kind(1.0d0)), dimension(n) :: y, func
  end function func
end interface

procedure(func), pointer :: pf => null()

contains

subroutine radau5_integration(f, tol, neq, n, y, tini, tend, info)
  implicit none
  interface
    function f(n, y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f
    end function f
  end interface
  real(kind=dp) :: tol
  integer :: neq
  integer :: n
  real(kind=dp), dimension(n) :: y
  real(kind=dp) :: tini, tend
  logical info
  ! local variables
  real(kind=dp) :: dt, t
  real(kind=dp) :: atol, rtol
  integer :: itol
  integer :: ijac, mljac, mujac
  integer :: imas, mlmas, mumas
  integer :: iout
  integer :: lwork
  real(kind=dp), allocatable, dimension(:) :: work   
  integer :: liwork
  integer, allocatable, dimension(:) :: iwork
  real(kind=dp) :: rpar
  integer :: ipar
  integer :: idid
  integer :: i, j
 
  ! initialisation
  pf => f

  ! initial time 
  t = tini

  ! initial step size  
  dt = 0.d0

  ! required tolerance
  rtol = tol
  atol = rtol
  itol = 0

  ! jacobian is computed internally by finite differences
  ijac = 0

  ! jacobian bandwidth
  mljac = neq
  mujac = neq

  ! mass-matrix is supposed to be the identity matrix
  ! subroutine mas is never called
  imas = 0
  mlmas = 0
  mumas = 0

  ! iout=0 : subroutine solout is never called
  ! iout=1 : subroutine solout is available for output
  iout = 0
  if (info) then
    iout = 1
  end if

  ! allocation of array of working space 
  ! case where the jacobian is full and the mass-matrix is the indentity 
  lwork = 4*n*n + 12*n + 20
  allocate(work(lwork))
  work(1:20) = 0.d0 
 
  ! allocation of array of integer working space
  liwork = 3*n + 20
  allocate(iwork(liwork))
  iwork(1:20) = 0 
  
  call radau5(n, f_radau, t, y, tend, dt, rtol, atol, itol,              &
              jac, ijac, mljac, mujac, mas, imas, mlmas, mumas,          &
              solout, iout, work, lwork, iwork, liwork, rpar, ipar, idid)

  if (info) then 
    print "(a, es12.3)", "  Tolerance : ", rtol
    print *, "  Number of computed steps :", iwork(16)
    print *, "  Number of accepted steps :", iwork(17)
    print *, "  Number of rejected steps :", iwork(18)
  endif  

  deallocate(work)
  deallocate(iwork)

end subroutine radau5_integration

subroutine f_radau(n, t, y, dy, rpar, ipar)
  implicit none
  ! parameter
  integer n
  real(kind=dp) :: t
  real(kind=dp), dimension(n) :: y, dy
  real(kind=dp) :: rpar
  integer :: ipar

  dy = pf(n, y)
end subroutine f_radau

subroutine solout(nr, xold, x, y, cont, lrc, n, rpar, ipar, irtrn)
  implicit none
  ! parameter
  integer nr
  real(kind=dp) :: xold, x
  integer :: n
  real(kind=dp), dimension(n) :: y
  integer :: lrc
  real(kind=dp), dimension(lrc) :: cont
  real(kind=dp) :: rpar
  integer :: ipar
  integer :: irtrn

  print *, x

end subroutine solout

subroutine jac()
end subroutine jac

subroutine mas()
end subroutine mas

end module mod_radau
